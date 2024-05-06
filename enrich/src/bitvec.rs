use uuid::Uuid;
use std::collections::HashMap;

// Representing indexes in a larger vector (background)
pub struct SparseBitVec{pub v:Vec<u32>}
impl SparseBitVec {
    pub fn new(background: &HashMap<Uuid, u32>, gene_set: &Vec<Uuid>) -> Self {
        SparseBitVec{v:gene_set.iter().filter_map(|gene_id| Some(*background.get(gene_id)?)).collect()}
    }
}

// Storing a multi-hot encoded vector for the gene set
#[derive(Eq, PartialEq, PartialOrd, Ord)]
pub struct DenseBitVec{pub v:Vec<u8>, pub n:usize}
impl DenseBitVec {
    pub fn new(background: &HashMap<Uuid, u32>, gene_set: &Vec<Uuid>) -> Self {
        let mut v: Vec<u8> = [0].repeat(background.len());
        let mut n = 0;
        for gene_id in gene_set {
            if let Some(gene_index) = background.get(&gene_id) {
                if v[*gene_index as usize] == 0 {
                    n += 1;
                    v[*gene_index as usize] = 1;
                }
            }
        }
        DenseBitVec{v, n}
    }
}

// this exploits the sparse vector to compute this overlap rather quickly
//  and in a way that is independent of the input gene set size
pub fn compute_overlap(a: &DenseBitVec, b: &SparseBitVec) -> u32 {
    b.v.iter().map(|i| a.v[*i as usize] as u32).sum()
}
