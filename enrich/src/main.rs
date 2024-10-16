mod async_rwlockhashmap;
mod fastfisher;
mod bitvec;

#[macro_use] extern crate rocket;
use async_lock::RwLock;
use futures::StreamExt;
use num::Integer;
use rocket::http::ContentType;
use std::future;
use std::io::Cursor;
use rocket::request::Request;
use rocket::response::{self, Response, Responder, stream::TextStream};
use rocket_db_pools::{Database, Connection};
use rocket_db_pools::sqlx::{self, Row};
use std::collections::HashMap;
use rayon::prelude::*;
use uuid::Uuid;
use adjustp::{adjust, Procedure};
use rocket::{State, response::status::Custom, http::Status};
use rocket::serde::{json::{json, Json, Value}, Serialize};
use std::sync::Arc;
use retainer::Cache;
use std::time::Instant;

use fastfisher::FastFisher;
use async_rwlockhashmap::RwLockHashMap;
use bitvec::{SparseBitVec,DenseBitVec,compute_overlap};

/**
 * Without this alternative allocator, very large chunks of memory do not get released back to the OS causing a large memory footprint over time.
 */
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[derive(Database)]
#[database("postgres")]
struct Postgres(sqlx::PgPool);

struct Bitmap<B: Integer + Copy + Into<usize>> {
    columns: HashMap<Uuid, B>,
    columns_str: Vec<String>,
    values: Vec<(Uuid, String, String, String, String, f64, SparseBitVec<B>)>,
}


impl<B: Integer + Copy + Into<usize>> Bitmap<B> {
    fn new() -> Self {
        Bitmap {
            columns: HashMap::new(),
            columns_str: Vec::new(),
            values: Vec::new(),
        }
    }
}

#[derive(Eq, PartialEq, PartialOrd, Ord)]
struct BackgroundQuery {
    background_id: Uuid,
    input_gene_set: DenseBitVec,
}

// This structure stores a persistent many-reader single-writer hashmap containing cached indexes for a given background id
struct PersistentState { 
    fisher: RwLock<FastFisher>,
    // NOTE: Bitmap<u16> limits the number of genes to 65K -- to support more than that, use u32/u64 at the cost of more memory
    bitmaps: RwLockHashMap<Uuid, Bitmap<u16>>,
    latest: RwLock<Option<Uuid>>,
    cache: Cache<Arc<BackgroundQuery>, Arc<Vec<PartialQueryResult>>>,
}

// The response data, containing the ids, and relevant metrics
#[derive(Clone)]
struct PartialQueryResult {
    index: usize,
    n_overlap: u32,
    odds_ratio: f64,
    pvalue: f64,
    adj_pvalue: f64,
    direction: String,
    n_gs_gene_id: u32,
    silhouette_score: f64,
}

#[derive(Serialize, Debug)]
struct QueryResult {
    gene_set_id: String,
    n_overlap: u32,
    odds_ratio: f64,
    pvalue: f64,
    adj_pvalue: f64,
}

struct QueryResponse {
    results: Vec<QueryResult>,
    content_range: (usize, usize, usize),
    terms: Vec<String>,
}

#[rocket::async_trait]
impl<'r> Responder<'r, 'static> for QueryResponse {
    fn respond_to(self, _: &'r Request<'_>) -> response::Result<'static> {
        let json = rocket::serde::json::serde_json::json!({
            "results": &self.results,
            "terms": &self.terms
        });
        let json_str = rocket::serde::json::serde_json::to_string(&json).unwrap();
        Response::build()
            .header(ContentType::JSON)
            .raw_header("Range-Unit", "items")
            .raw_header("Content-Range", format!("{}-{}/{}", self.content_range.0, self.content_range.1, self.content_range.2))
            .sized_body(json_str.len(), Cursor::new(json_str))
            .ok()
    }
}


// Ensure the specific background_id exists in state, resolving it if necessary
async fn ensure_index(db: &mut Connection<Postgres>, state: &State<PersistentState>, background_id: Uuid) -> Result<(), String> {
    if state.bitmaps.contains_key(&background_id).await {
        return Ok(())
    }

    println!("[{}] initializing", background_id);
    let start = Instant::now();
    {
        // this lets us write a new bitmap by only blocking the whole hashmap for a short period to register the new bitmap
        // after which we block the new empty bitmap for writing
        let mut bitmap = state.bitmaps.insert_write(background_id, Bitmap::new()).await;

        let background_info = sqlx::query("select id, species, (select jsonb_object_agg(g.id, g.symbol) from jsonb_each(gene_ids) bg(gene_id, nil) inner join app_public_v2.gene g on bg.gene_id::uuid = g.id) as genes from app_public_v2.background b where id = $1::uuid;")
            .bind(background_id.to_string())
            .fetch_one(&mut **db).await.map_err(|e| e.to_string())?;

        let background_genes: sqlx::types::Json<HashMap<String, String>> = background_info.try_get("genes").map_err(|e| e.to_string())?;
        let species: String = background_info.try_get("species").map_err(|e| e.to_string())?;
        let mut background_genes = background_genes.iter().map(|(id, symbol)| Ok((Uuid::parse_str(id).map_err(|e| e.to_string())?, symbol.clone()))).collect::<Result<Vec<_>, String>>()?;
        background_genes.sort_unstable();
        {
            let mut fisher = state.fisher.write().await;
            fisher.extend_to(background_genes.len()*4);
        };
        bitmap.columns.reserve(background_genes.len());
        bitmap.columns_str.reserve(background_genes.len());
        for (i, (gene_id, gene)) in background_genes.into_iter().enumerate() {
            bitmap.columns.insert(gene_id, i as u16);
            bitmap.columns_str.push(gene);
        }

        // compute the index in memory
        sqlx::query(format!("select gs.id, gs.term, gs.gene_ids, gp.title, gp.silhouette_score, gp.gse_attrs, et.sig_terms from app_public_v2.gene_set gs join app_public_v2.gene_set_pmid gp on gs.id = gp.id left join app_public_v2.enrichr_terms et on gp.term = et.sig where gs.species = '{}';", species).as_str())
            .fetch(&mut **db)
            .for_each(|row| {
                let row = row.unwrap();
                let gene_set_id: uuid::Uuid = row.try_get("id").unwrap();
                let term: String = row.try_get("term").unwrap();
                let title: String = row.try_get("title").unwrap();
                let attrs: Result<Option<String>, _> = row.try_get("gse_attrs");
                let attrs = attrs.unwrap_or_else(|_| None).unwrap_or_default();
                let sig_terms: Option<Vec<String>> = row.try_get("sig_terms").unwrap_or_default();
                let sig_terms_str = sig_terms.map(|terms| terms.join(" ")).unwrap_or_else(|| String::new());
                let silhouette_score: f64 = row.try_get("silhouette_score").unwrap_or_default();
                let gene_ids: sqlx::types::Json<HashMap<String, sqlx::types::JsonValue>> = row.try_get("gene_ids").unwrap();
                let gene_ids = gene_ids.keys().map(|gene_id| Uuid::parse_str(gene_id).unwrap()).collect::<Vec<Uuid>>();
                let bitset = SparseBitVec::new(&bitmap.columns, &gene_ids);
                bitmap.values.push((gene_set_id, term, title, attrs, sig_terms_str, silhouette_score, bitset));
                future::ready(())
            })
            .await;
    }
    let duration = start.elapsed();
    println!("[{}] initialized in {:?}", background_id, duration);
    {
        let mut latest = state.latest.write().await;
        latest.replace(background_id);
    }
    Ok(())
}

#[get("/<background_id>")]
async fn ensure(
    mut db: Connection<Postgres>,
    state: &State<PersistentState>,
    background_id: &str,
) -> Result<Value, Custom<String>> {
    let background_id = Uuid::parse_str(background_id).map_err(|e| Custom(Status::BadRequest, e.to_string()))?;
    ensure_index(&mut db, &state, background_id).await.map_err(|e| Custom(Status::InternalServerError, e.to_string()))?;
    let bitmap = state.bitmaps.get_read(&background_id).await.ok_or(Custom(Status::NotFound, String::from("Can't find background")))?;
    Ok(json!({
        "columns": bitmap.columns.len(),
        "index": bitmap.values.len(),
    }))
}

/**
 * This is a helper for building a GMT file on the fly, it's much cheaper to do this here
 *  than fetch it from the database, it's also nice since we won't need to save raw files.
 */
#[get("/<background_id>/gmt")]
async fn get_gmt(
    mut db: Connection<Postgres>,
    state: &State<PersistentState>,
    background_id: String,
) -> Result<TextStream![String + '_], Custom<String>> {
    let background_id = {
        if background_id == "latest" {
            let latest = state.latest.read().await;
            latest.clone().ok_or(Custom(Status::NotFound, String::from("Nothing loaded")))?
        } else {
            Uuid::parse_str(&background_id).map_err(|e| Custom(Status::BadRequest, e.to_string()))?
        }
    };
    ensure_index(&mut db, &state, background_id).await.map_err(|e| Custom(Status::InternalServerError, e.to_string()))?;
    let bitmap = state.bitmaps.get_read(&background_id).await.ok_or(Custom(Status::InternalServerError, String::from("Can't find background")))?;
    Ok(TextStream! {
        for (_row_id, row_str, _title, _attrs, _sig_terms, _silhouette_score, gene_set) in bitmap.values.iter() {
            let mut line = String::new();
            line.push_str(row_str);
            line.push_str("\t");
            for col_ind in gene_set.v.iter() {
                line.push_str("\t");
                line.push_str(&bitmap.columns_str[*col_ind as usize]);
            }
            line.push_str("\n");
            yield line
        }
    })
}

#[delete("/<background_id>")]
async fn delete(
    state: &State<PersistentState>,
    background_id: &str,
) -> Result<(), Custom<String>> {
    let background_id = {
        if background_id == "latest" {
            let latest = state.latest.read().await;
            latest.clone().ok_or(Custom(Status::NotFound, String::from("Nothing loaded")))?
        } else {
            Uuid::parse_str(&background_id).map_err(|e| Custom(Status::BadRequest, e.to_string()))?
        }
    };
    if !state.bitmaps.contains_key(&background_id).await {
        return Err(Custom(Status::NotFound, String::from("Not Found")));
    }
    if state.bitmaps.remove(&background_id).await {
        println!("[{}] deleted", background_id);
    }
    Ok(())
}

// query a specific background_id, providing the bitset vector as input
//  the result are the gene_set_ids & relevant metrics
// this can be pretty fast since the index is saved in memory and the overlaps can be computed in parallel
#[post("/<background_id>?<filter_term>&<overlap_ge>&<pvalue_le>&<adj_pvalue_le>&<offset>&<limit>&<score_filter>&<sort_by>&<sort_by_dir>", data = "<input_gene_set>")]
async fn query(
    mut db: Connection<Postgres>,
    state: &State<PersistentState>,
    input_gene_set: Json<Vec<String>>,
    background_id: &str,
    filter_term: Option<String>,
    overlap_ge: Option<u32>,
    pvalue_le: Option<f64>,
    adj_pvalue_le: Option<f64>,
    offset: Option<usize>,
    limit: Option<usize>,
    score_filter: Option<f64>,
    sort_by: Option<String>,
    sort_by_dir: Option<String>
) -> Result<QueryResponse, Custom<String>> {
    let background_id = {
        if background_id == "latest" {
            let latest = state.latest.read().await;
            latest.clone().ok_or(Custom(Status::NotFound, String::from("Nothing loaded")))?
        } else {
            Uuid::parse_str(&background_id).map_err(|e| Custom(Status::BadRequest, e.to_string()))?
        }
    };
    ensure_index(&mut db, &state, background_id).await.map_err(|e| Custom(Status::InternalServerError, e.to_string()))?;
    let start = Instant::now();
    let input_gene_set = input_gene_set.0.into_iter().map(|gene| Uuid::parse_str(&gene)).collect::<Result<Vec<_>, _>>().map_err(|e| Custom(Status::BadRequest, e.to_string()))?;
    let bitmap = state.bitmaps.get_read(&background_id).await.ok_or(Custom(Status::NotFound, String::from("Can't find background")))?;
    let input_gene_set = DenseBitVec::new(&bitmap.columns, &input_gene_set);
    let filter_term = filter_term.and_then(|filter_term| Some(filter_term.to_lowercase()));
    let overlap_ge = overlap_ge.unwrap_or(1);
    let pvalue_le =  pvalue_le.unwrap_or(1.0);
    let adj_pvalue_le =  adj_pvalue_le.unwrap_or(1.0);
    let score_filter = score_filter.unwrap_or(-1.0);
    let sort_by = sort_by.unwrap_or("pvalue".to_string());
    let sort_by_dir = sort_by_dir.unwrap_or("asc".to_string());
    let background_query = Arc::new(BackgroundQuery { background_id, input_gene_set });
    let results = {
        let results = state.cache.get(&background_query).await;
        if let Some(results) = results {
            results.value().clone()
        } else {
            // parallel overlap computation
            let n_background = bitmap.columns.len() as u32;
            let n_user_gene_id = background_query.input_gene_set.n as u32;
            let fisher = state.fisher.read().await;
            let mut results: Vec<_> = bitmap.values.par_iter()
                .enumerate()
                .filter_map(|(index, (_row_id, row_str, _title_str, _attrs, _sig_terms, silhouette_score, gene_set))| {
                    let n_overlap = compute_overlap(&background_query.input_gene_set, &gene_set) as u32;
                    if n_overlap < overlap_ge {
                        return None
                    }
                    let n_gs_gene_id = gene_set.v.len() as u32;
                    let a = n_overlap;
                    let b = n_user_gene_id - a;
                    let c = n_gs_gene_id - a;
                    let d = n_background - b - c + a;
                    let pvalue = fisher.get_p_value(a as usize, b as usize, c as usize, d as usize);
                    if pvalue > pvalue_le {
                        return None
                    }
                    let direction = row_str.split_whitespace().last().unwrap_or("").to_string();
                    let odds_ratio = ((n_overlap as f64) / (n_user_gene_id as f64)) / ((n_gs_gene_id as f64) / (n_background as f64));
                    Some(PartialQueryResult { index, n_overlap, odds_ratio, pvalue ,adj_pvalue: 1.0, direction, n_gs_gene_id, silhouette_score: *silhouette_score })
                })
                .collect();
            // extract pvalues from results and compute adj_pvalues
            let mut pvalues = vec![1.0; bitmap.values.len()];
            for result in &results {
                pvalues[result.index] = result.pvalue;
            }
            // add adj_pvalues to results
            let adj_pvalues = adjust(&pvalues, Procedure::BenjaminiHochberg);
            results.retain_mut(|result| {
                if let Some(adj_pvalue) = adj_pvalues.get(result.index) {
                    result.adj_pvalue = *adj_pvalue;
                }
                result.adj_pvalue <= adj_pvalue_le
            });

            let results = Arc::new(results);
            state.cache.insert(background_query, results.clone(), 30000).await;
            let duration = start.elapsed();
            println!("[{}] {} genes enriched in {:?} sorted by {} {}", background_id, n_user_gene_id, duration, sort_by, sort_by_dir);
            results
        }
    };
    // Assuming `results` is of type Arc<Vec<PartialQueryResult>>
    // Clone the data out of the Arc for sorting
    let mut results_vec: Vec<PartialQueryResult> = (*results).clone();

    // Perform the sorting on the cloned data
    results_vec.sort_unstable_by(|a, b| {
        let order = match sort_by.as_str() {
            "pvalue" => a.pvalue.partial_cmp(&b.pvalue).unwrap_or(std::cmp::Ordering::Equal),
            "adj_pvalue" => a.adj_pvalue.partial_cmp(&b.adj_pvalue).unwrap_or(std::cmp::Ordering::Equal),
            "odds_ratio" => a.odds_ratio.partial_cmp(&b.odds_ratio).unwrap_or(std::cmp::Ordering::Equal),
            "n_overlap" => a.n_overlap.partial_cmp(&b.n_overlap).unwrap_or(std::cmp::Ordering::Equal),
            "size" => a.n_gs_gene_id.partial_cmp(&b.n_gs_gene_id).unwrap_or(std::cmp::Ordering::Equal),
            "silhouette_score" => a.silhouette_score.partial_cmp(&b.silhouette_score).unwrap_or(std::cmp::Ordering::Equal),
            "direction" => a.direction.partial_cmp(&b.direction).unwrap_or(std::cmp::Ordering::Equal),
            _ => a.pvalue.partial_cmp(&b.pvalue).unwrap_or(std::cmp::Ordering::Equal),
        };
        // Adjust the sorting direction based on sort_by_dir
        match sort_by_dir.as_str() {
            "desc" => order.reverse(),
            _ => order,
        }
    });

    // Wrap the sorted data back into a new Arc
    let sorted_results = Arc::new(results_vec);
    let mut all_enriched_terms: Vec<String> = Vec::new();
    let mut results: Vec<_> = sorted_results
        .iter()
        .filter_map(|result| {
            let (gene_set_id, gene_set_term, title, attrs, sig_terms, silhouette_score, _gene_set) = bitmap.values.get(result.index)?;
            if let Some(filter_term) = &filter_term {
                if !title.to_lowercase().contains(filter_term) && !gene_set_term.to_lowercase().contains(filter_term) && !attrs.contains(filter_term) && !sig_terms.to_lowercase().contains(filter_term) { 
                    return None 
                }

            } 
            if *silhouette_score < score_filter {
                return None
            }
            all_enriched_terms.push(gene_set_term.clone());
            Some(QueryResult {
                gene_set_id: gene_set_id.to_string(),
                n_overlap: result.n_overlap,
                odds_ratio: result.odds_ratio,
                pvalue: result.pvalue,
                adj_pvalue: result.adj_pvalue,
            })
        })
        .collect();
    let range_total = results.len();

    let (range_start, range_end) = match (offset.unwrap_or(0), limit) {
        (0, None) => (0, range_total),
        (offset, None) => {
            if offset < results.len() {
                results.drain(..offset);
            };
            (offset, range_total)
        },
        (offset, Some(limit)) => {
            if offset < results.len() {
                results.drain(..offset);
                if limit < results.len() {
                    results.drain(limit..);
                }
            };
            (offset, offset + results.len())
        },
    };
    Ok(QueryResponse {
        results,
        content_range: (range_start, range_end, range_total),
        terms: all_enriched_terms,
    })
}

#[launch]
fn rocket() -> _ {
    rocket::build()
        .manage(PersistentState {
            fisher: RwLock::new(FastFisher::new()),
            bitmaps: RwLockHashMap::new(),
            latest: RwLock::new(None),
            cache: Cache::new(),
        })
        .attach(Postgres::init())
        .mount("/", routes![ensure, get_gmt, query, delete])
}
