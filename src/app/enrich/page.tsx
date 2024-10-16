"use client";
import React from "react";
import {
  FetchUserGeneSetQuery,
  useEnrichmentQueryQuery,
  useFetchUserGeneSetQuery,
  useOverlapQueryQuery,
  useViewGeneSetQuery,
  useGetBackgroundsQuery,
  EnrichedTermResult,
  useTermEnrichmentQuery,
  useEnrichrTermEnrichmentQuery,
  EnrichrResult,
} from "@/graphql";
import determineSpecies from "@/utils/determineSpecies";
import HypothesisModal from "@/components/hypothesisModal";
import ensureArray from "@/utils/ensureArray";
import Loading from "@/components/loading";
import Pagination from "@/components/pagination";
import useQsState from "@/utils/useQsState";
import Stats from "../stats";
import Image from "next/image";
import GeneSetModal from "@/components/geneSetModal";
import TermVis from "@/components/termVis";
import EnrichrTermVis from "@/components/enrichrTermVis";
import SamplesModal from "@/components/samplesModal";
import EnrichrModal from "@/components/enrichrModal";
import partition from "@/utils/partition";
import { FaInfo, FaSearch } from "react-icons/fa";
import { TiDeleteOutline } from "react-icons/ti";
import { FaFilter } from "react-icons/fa";
import { MdOutlineFileDownload } from "react-icons/md";
import { RiAiGenerate } from "react-icons/ri";
import { FaSort } from "react-icons/fa";
import classNames from "classnames";

const pageSize = 8;

type GeneSetModalT =
  | {
      type: "UserGeneSet";
      description: string;
      genes: string[];
    }
  | {
      type: "GeneSetOverlap";
      id: string;
      description: string;
      genes: string[];
      gseId?: string;
    }
  | {
      type: "GeneSet";
      id: string;
      description: string;
    }
  | undefined;

function EnrichmentResults({
  userGeneSet,
  setModalGeneSet,
  setModalHypothesis,
  setModalSamples,
  setModalCondition,
  setModalEnrichr
}: {
  userGeneSet?: FetchUserGeneSetQuery;
  setModalGeneSet: React.Dispatch<React.SetStateAction<GeneSetModalT>>;
  setModalHypothesis: React.Dispatch<React.SetStateAction<GeneSetModalT>>;
  setModalSamples: React.Dispatch<React.SetStateAction<string[] | undefined>>;
  setModalCondition: React.Dispatch<React.SetStateAction<string | undefined>>;
  setModalEnrichr: React.Dispatch<React.SetStateAction<GeneSetModalT>>;
}) {
  const genes = React.useMemo(
    () =>
      ensureArray(userGeneSet?.userGeneSet?.genes).filter(
        (gene): gene is string => !!gene
      ),
    [userGeneSet]
  );
  const species = React.useMemo(
    () => determineSpecies(genes[0] || ""),
    [genes]
  );
  const [tab, setTab] = React.useState(1);
  const [sourceType, setSourceType] = React.useState("A");

  const { data: backgrounds } = useGetBackgroundsQuery();
  var backgroundIds: Record<string, string> = {};
  backgrounds?.backgrounds?.nodes?.forEach((background) => {
    backgroundIds[background?.species ?? ""] = background?.id ?? "";
  });
  const [queryString, setQueryString] = useQsState({ page: "1", q: "", sortBy: "pvalue", sortByDir: "asc"});
  const [rawTerm, setRawTerm] = React.useState("");
  const [enrichedTerms, setEnrichedTerms] = React.useState<(string | null)[]>();
  const [topEnrichedSigs, setTopEnrichedSigs] =
    React.useState<(string | null)[]>();
  const [filterScore, setFilterScore] = React.useState(-1);
  const [filterScoreSilder, setFilterScoreSilder] = React.useState(-1);
  const { page, term, sortBy, sortByDir } = React.useMemo(
    () => ({
      page: queryString.page ? + queryString.page : 1,
      term: queryString.q ?? "",
      sortBy: queryString.sortBy ?? "pvalue",
      sortByDir: queryString.sortByDir ?? "asc",
    }),
    [queryString]
  );
  const { data: enrichmentResults, refetch, error } = useEnrichmentQueryQuery({
    skip: genes.length === 0,
    variables: {
      genes,
      filterTerm: term,
      offset: (page - 1) * pageSize,
      first: pageSize,
      id: backgroundIds[species],
      filterScoreLe: filterScore,
      sortBy: sortBy,
      sortByDir: sortByDir,
    },
  });

  console.log(enrichmentResults)

  React.useEffect(() => {
    setRawTerm(term);
  }, [term]);

  React.useEffect(() => {
    if (enrichmentResults?.background?.enrich?.enrichedTerms) {
      setEnrichedTerms(enrichmentResults?.background?.enrich?.enrichedTerms);
    }
  }, [enrichmentResults, setEnrichedTerms]);

  React.useEffect(() => {
    if (enrichmentResults?.background?.enrich?.topEnrichedSigs) {
      setTopEnrichedSigs(
        enrichmentResults?.background?.enrich?.topEnrichedSigs
      );
    }
  }, [enrichmentResults, setTopEnrichedSigs]);

  const { data: termEnrichmentResults } = useTermEnrichmentQuery({
    variables: {
      enrichedTerms: enrichedTerms,
      sourceType: sourceType,
      species: species,
    },
  });

  const { data: EnrichrTermEnrichmentResults } = useEnrichrTermEnrichmentQuery({
    variables: {
      enrichedSigs: topEnrichedSigs,
      species: species,
    },
  });

  function setSort(sortByCol: string) {
    if (sortByCol === sortBy) {
      setQueryString({page: "1", q: term, sortByDir: sortByDir === "asc" ? "desc" : "asc", sortBy: sortByCol});
      refetch({sortByDir: sortByDir === "asc" ? "desc" : "asc" ,sortBy: sortByCol})
    } else {
      setQueryString({page: "1", q: term, sortByDir: sortByDir, sortBy: sortByCol});
      refetch({sortByDir:sortByDir, sortBy: sortByCol})
    } 
  }

  if (error) {
    return <div className="m-5">No valid genes identified.</div>
  }

  return (
    <div className="flex flex-col gap-2 my-2">
      <ul
        className="relative flex flex-wrap p-1 list-none rounded-lg bg-inherit gap-1"
        data-tabs="tabs"
        role="list"
      >
        <li
          className={classNames(
            "z-30 flex-auto text-center p-2 cursor-pointer border rounded-md",
            {
              "font-bold text-white bg-slate-700 bg-opacity-50 rounded-lg":
                tab === 1,
            }
          )}
        >
          <a
            className="z-30 flex items-center justify-center w-full px-0 py-1 mb-0 transition-all ease-in-out rounded-lg cursor-pointerbg-inherit"
            onClick={(evt) => {
              evt.preventDefault();
              setTab(1);
            }}
          >
            <span className="ml-1">Matching Gene Sets</span>
          </a>
        </li>
        <li
          className={classNames(
            "z-30 flex-auto text-center p-2 cursor-pointer border rounded-md",
            {
              "font-bold text-white bg-slate-700 bg-opacity-50 rounded-lg":
                tab === 2,
            }
          )}
        >
          <a
            className="z-30 flex items-center justify-center w-full px-0 py-1 mb-0 transition-all ease-in-out rounded-lg cursor-pointerbg-inherit"
            onClick={(evt) => {
              evt.preventDefault();
              setTab(2);
            }}
          >
            <span className="ml-1">Common Terms in Matching Gene Sets</span>
          </a>
        </li>
        <li
          className={classNames(
            "z-30 flex-auto text-center p-2 cursor-pointer border rounded-md",
            {
              "font-bold text-white bg-slate-700 bg-opacity-50 rounded-lg":
                tab === 3,
            }
          )}
        >
          <a
            className="z-30 flex items-center justify-center w-full px-0 py-1 mb-0 transition-all ease-in-out rounded-lg cursor-pointerbg-inherit"
            onClick={(evt) => {
              evt.preventDefault();
              setTab(3);
            }}
          >
            <span className="ml-1">Enrichr Terms</span>
          </a>
        </li>
      </ul>

      <h2 className="text-md font-bold">
        {!enrichmentResults?.background?.enrich ? (
          <>
            Rummaging through{" "}
            {species == "human" ? (
              <>
                <Stats show_human_gene_sets /> gene sets
              </>
            ) : (
              <Stats show_mouse_gene_sets />
            )}{" "}
          </>
        ) : (
          <>
            After rummaging through{" "}
            {species == "human" ? (
              <Stats show_human_gene_sets />
            ) : (
              <Stats show_mouse_gene_sets />
            )}{" "}
            gene sets. RummaGEO{" "}
            <Image
              className="inline-block rounded"
              src="/images/rummageo_logo.png"
              width={50}
              height={100}
              alt="Rummageo"
            ></Image>{" "}
            found{" "}
            {Intl.NumberFormat("en-US", {}).format(
              enrichmentResults?.background?.enrich?.totalCount || 0
            )}{" "}
            statistically significant matches.
          </>
        )}
      </h2>
      {!enrichmentResults?.background?.enrich ? <Loading /> : null}
      {tab == 1 && enrichmentResults ? (
        <>
          <div className="flex flex-wrap items-center justify-end gap-3">
            <p className="text-sm">
              Minimum Silhouette Score: <b>{filterScoreSilder}</b>
            </p>
            <input
              id="default-range"
              type="range"
              value={filterScoreSilder}
              onChange={(evt) => setFilterScoreSilder(Number(evt.target.value))}
              max={1}
              min={-1}
              step={0.01}
              className="w-1/8 h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer dark:bg-gray-700"
            ></input>
            <div className="tooltip" data-tip="Filter results">
              <button
                className="btn bg-transparent"
                onClick={() => {
                  setQueryString({ page: "1", q: "", sortBy: sortBy, sortByDir: sortByDir})
                  setFilterScore(filterScoreSilder)
                }
                  
                }
              >
                <FaFilter />
              </button>
            </div>
            <form
              className="join flex flex-row place-content-end place-items-center"
              onSubmit={(evt) => {
                evt.preventDefault();
                setQueryString({ page: "1", q: rawTerm, sortBy: sortBy, sortByDir: sortByDir});
              }}
            >
              <input
                type="text"
                className="input input-bordered bg-transparent join-item"
                value={rawTerm}
                onChange={(evt) => {
                  setRawTerm(evt.currentTarget.value);
                }}
              />
              <div className="tooltip" data-tip="Search results">
                <button
                  type="submit"
                  className="btn join-item bg-transparent ml-2"
                >
                  <FaSearch />
                </button>
              </div>
              <div className="tooltip" data-tip="Clear search">
                <button
                  type="reset"
                  className="btn join-item bg-transparent"
                  onClick={(evt) => {
                    setQueryString({ page: "1", q: "",  sortBy: sortBy, sortByDir: sortByDir});
                  }}
                >
                  <TiDeleteOutline />
                </button>
              </div>
              <a
                href={`/enrich/download?dataset=${queryString.dataset}&q=${queryString.q}`}
                download="results.tsv"
              >
                <div className="tooltip" data-tip="Download results">
                  <button
                    type="button"
                    className="btn join-item font-bold text-2xl pb-1 bg-transparent"
                  >
                    <MdOutlineFileDownload />
                  </button>
                </div>
              </a>
            </form>
            <div
              className="tooltip z-40"
              data-tip="
            Submit a description of your gene set to generate a hypothesis utlizing the signature's associated study 
            abstract or summary and highly enriched terms from Enrichr libraries: WikiPathway 2023 Human GWAS_Catalog 2023, 
            GO Biological Process 2023, and MGI Mammalian Phenotype Level 4 2021. "
            >
              <div className="bg-transparent text-xs h-fit p-2 flex text-center border rounded-lg mx-auto">
                <p className="font-light">Hypothesis Generation</p> <FaInfo className="ml-1 mt-0.5"/>
              </div>
            </div>
          </div>
          <div className="border-2 rounded-lg p-5">
            <div className="overflow-x-auto">
              <table className="table table-xs">
                <thead>
                  <tr>
                    <th>GEO Series</th>
                    <th className="hidden  2xl:table-cell">PMID</th>
                    <th>Title</th>
                    <th>Condition 1</th>
                    <th>Condition 2</th>
                    <th><span className="flex align-text-top cursor-pointer" onClick={() => setSort("direction")}>Direction <FaSort size={15} /></span></th>
                    <th className="hidden  2xl:table-cell">Platform</th>
                    <th className="hidden  2xl:table-cell">Date</th>
                    <th><span className="flex align-text-top cursor-pointer" onClick={() => setSort("size")}>Gene Set Size <FaSort size={15}  /></span></th>
                    <th ><span className="flex align-text-top cursor-pointer" onClick={() => setSort("n_overlap")}>Overlap<FaSort size={15}/></span></th>
                    <th className="hidden  2xl:table-cell"><span className="flex align-text-top cursor-pointer" onClick={() => setSort("odds_ratio")}>Odds<FaSort size={15}/></span></th>
                    <th className="hidden  2xl:table-cell"><span className="flex align-text-top cursor-pointer" onClick={() => setSort("pvalue")}>PValue<FaSort size={15}/></span></th>
                    <th><span className="flex align-text-top cursor-pointer" onClick={() => setSort("pvalue")}>Adj. PValue<FaSort size={15}/></span></th>
                    <th className="hidden  2xl:table-cell"><span className="flex align-text-top cursor-pointer" onClick={() => setSort("silhouette_score")}>Silhouette Score<FaSort size={15}/></span></th>
                    <th>Hypothesis</th>
                    <th>Enrichr Terms</th>
                  </tr>
                </thead>
                <tbody>
                  {enrichmentResults?.background?.enrich?.nodes?.map(
                    (enrichmentResult, j) => {
                      if (!enrichmentResult?.geneSet) return null;
                      const [gse, cond1, _, cond2, __, dir] = partition(
                        enrichmentResult?.geneSet?.term
                      );

                      const term = enrichmentResult?.geneSet?.term;
                      var pmid =
                        enrichmentResult?.geneSet?.geneSetPmidsById?.nodes[0]
                          ?.pmid ?? null;
                      if (pmid?.includes(",")) {
                        pmid = JSON.parse(pmid.replace(/'/g, '"')).join(",");
                      }
                      var platform =
                        enrichmentResult?.geneSet?.geneSetPmidsById?.nodes[0]
                          ?.platform ?? "";
                      if (platform?.includes(",")) {
                        platform = JSON.parse(platform.replace(/'/g, '"')).join(
                          ","
                        );
                      }

                      var node;

                      if (
                        enrichmentResult.geneSet?.geneSetPmidsById?.nodes
                          .length > 1
                      ) {
                        const sampleKeys = Object.keys(
                          enrichmentResult.geneSet?.geneSetPmidsById?.nodes[0]
                            ?.sampleGroups.samples
                        );

                        if (
                          sampleKeys.includes(cond1) &&
                          sampleKeys.includes(cond2)
                        ) {
                          node =
                            enrichmentResult.geneSet?.geneSetPmidsById
                              ?.nodes[0];
                        } else {
                          node =
                            enrichmentResult.geneSet?.geneSetPmidsById
                              ?.nodes[1];
                        }
                      } else {
                        node =
                          enrichmentResult.geneSet?.geneSetPmidsById?.nodes[0];
                      }

                      const cond1Title =
                        node?.sampleGroups?.titles[cond1] ?? "";
                      const cond2Title =
                        node?.sampleGroups?.titles[cond2] ?? "";

                      const cond1Samples =
                        node?.sampleGroups?.samples[cond1] ?? "";
                      const cond2Samples =
                        node?.sampleGroups?.samples[cond2] ?? "";

                      return (
                        <tr key={j} className="text-center">
                          <th>
                            {gse.includes(",") ? (
                              <>
                                {gse.split(",").map((g, i) => {
                                  return (
                                    <>
                                      <a
                                        key={i}
                                        className="underline cursor-pointer"
                                        href={`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${g}`}
                                        target="_blank"
                                        rel="noreferrer"
                                      >
                                        {g}
                                      </a>
                                      {i != gse.split(",").length - 1 ? (
                                        <>,</>
                                      ) : (
                                        <></>
                                      )}{" "}
                                    </>
                                  );
                                })}
                              </>
                            ) : (
                              <a
                                className="underline cursor-pointer"
                                href={`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gse}`}
                                target="_blank"
                                rel="noreferrer"
                              >
                                {gse}
                              </a>
                            )}
                          </th>
                          <th className="hidden  2xl:table-cell">
                            {pmid ? (
                              pmid.includes(",") ? (
                                <>
                                  {pmid.split(",").map((p, i) => {
                                    return (
                                      <>
                                        <a
                                          key={i}
                                          className="underline cursor-pointer"
                                          href={`https://pubmed.ncbi.nlm.nih.gov/${p}/`}
                                          target="_blank"
                                          rel="noreferrer"
                                        >
                                          {p}
                                        </a>
                                        {pmid ? (
                                          i != pmid?.split(",")?.length - 1 ? (
                                            <>,</>
                                          ) : (
                                            <></>
                                          )
                                        ) : (
                                          <></>
                                        )}{" "}
                                      </>
                                    );
                                  })}{" "}
                                </>
                              ) : (
                                <a
                                  className="underline cursor-pointer"
                                  href={`https://pubmed.ncbi.nlm.nih.gov/${pmid}/`}
                                  target="_blank"
                                  rel="noreferrer"
                                >
                                  {pmid}
                                </a>
                              )
                            ) : (
                              <>N/A</>
                            )}
                          </th>
                          <td className="text-left">
                            {enrichmentResult?.geneSet?.geneSetPmidsById
                              ?.nodes[0]?.title ?? ""}
                          </td>
                          <td className="text-left">
                            <label
                              htmlFor="geneSetModal"
                              className="prose underline cursor-pointer"
                              onClick={(evt) => {
                                setModalSamples(cond1Samples);
                                setModalCondition(cond1Title);
                              }}
                            >
                              {cond1Title}
                            </label>
                          </td>
                          <td className="text-left">
                            <label
                              htmlFor="geneSetModal"
                              className="prose underline cursor-pointer"
                              onClick={(evt) => {
                                setModalSamples(cond2Samples);
                                setModalCondition(cond2Title);
                              }}
                            >
                              {cond2Title}
                            </label>
                          </td>
                          <td>
                            {dir === "up"
                              ? "Up"
                              : dir === "dn"
                              ? "Down"
                              : "Up/Down"}
                          </td>
                          <td className="hidden  2xl:table-cell">
                            {platform ? (
                              platform.includes(",") ? (
                                <>
                                  {platform.split(",").map((p, i) => {
                                    return (
                                      <>
                                        <a
                                          key={i}
                                          className="underline cursor-pointer"
                                          href={`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${p}`}
                                          target="_blank"
                                          rel="noreferrer"
                                        >
                                          {p}
                                        </a>
                                        {i != platform.split(",").length - 1 ? (
                                          <>,</>
                                        ) : (
                                          <></>
                                        )}{" "}
                                      </>
                                    );
                                  })}{" "}
                                </>
                              ) : (
                                <a
                                  className="underline cursor-pointer"
                                  href={`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${platform}`}
                                  target="_blank"
                                  rel="noreferrer"
                                >
                                  {platform}
                                </a>
                              )
                            ) : (
                              <>N/A</>
                            )}
                          </td>
                          <td className="hidden  2xl:table-cell">
                            {enrichmentResult?.geneSet?.geneSetPmidsById
                              ?.nodes[0]?.publishedDate ?? ""}
                          </td>
                          <td className="whitespace-nowrap text-underline cursor-pointer">
                            <label
                              htmlFor="geneSetModal"
                              className="prose underline cursor-pointer"
                              onClick={(evt) => {
                                setModalGeneSet({
                                  type: "GeneSet",
                                  id: enrichmentResult?.geneSet?.id,
                                  description:
                                    enrichmentResult?.geneSet?.term ?? "",
                                });
                              }}
                            >
                              {enrichmentResult?.geneSet?.nGeneIds}
                            </label>
                          </td>
                          <td className="whitespace-nowrap text-underline cursor-pointer">
                            <label
                              htmlFor="geneSetModal"
                              className="prose underline cursor-pointer"
                              onClick={(evt) => {
                                setModalGeneSet({
                                  type: "GeneSetOverlap",
                                  id: enrichmentResult?.geneSet?.id,
                                  description:
                                    enrichmentResult?.geneSet?.term ?? "",
                                  genes,
                                });
                              }}
                            >
                              {enrichmentResult?.nOverlap}
                            </label>
                          </td>
                          <td className="whitespace-nowrap hidden 2xl:table-cell">
                            {enrichmentResult?.oddsRatio?.toPrecision(3)}
                          </td>
                          <td className="whitespace-nowrap hidden 2xl:table-cell">
                            {enrichmentResult?.pvalue?.toPrecision(3)}
                          </td>
                          <td className="whitespace-nowrap">
                            {enrichmentResult?.adjPvalue?.toExponential(2)}
                          </td>
                          <td className="whitespace-nowrap hidden 2xl:table-cell">
                            {enrichmentResult?.geneSet?.geneSetPmidsById?.nodes[0]?.silhouetteScore?.toPrecision(
                              2
                            ) ?? ""}
                          </td>
                          <td>
                            <div
                              className="tooltip tooltip-left"
                              data-tip="Generate GPT-4 Hypothesis"
                            >
                              <button
                                className="btn btn-sm"
                                onClick={(evt) => {
                                  setModalHypothesis({
                                    type: "GeneSetOverlap",
                                    id: enrichmentResult?.geneSet?.id,
                                    description:
                                      `${gse}: ${cond1Title} vs. ${cond1Title} ${dir}` ??
                                      "",
                                    genes,
                                    gseId:
                                      enrichmentResult?.geneSet
                                        ?.geneSetPmidsById.nodes[0]?.gseId ??
                                      "",
                                  });
                                }}
                              >
                                <RiAiGenerate />
                              </button>
                            </div>
                          </td>
                          <td>
                            <div
                              className="tooltip tooltip-left"
                              data-tip="View Enrichr Terms"
                            >
                              <button
                                className="btn btn-sm"
                                onClick={(evt) => {
                                  setModalEnrichr({
                                    type: "GeneSet",
                                    id: enrichmentResult?.geneSet?.id,
                                    description: term ?? "",
                                  });
                                }}
                              >
                                <Image
                                  src="/images/enrichr-icon.png"
                                  alt="Enrichr"
                                  width={20}
                                  height={20}
                                  className=""/>
                              </button>
                            </div>
                          </td>
                        </tr>
                      );
                    }
                  )}
                </tbody>
              </table>
            </div>
          </div>
          {enrichmentResults?.background?.enrich ? (
            <div className="w-full flex flex-col items-center">
              <Pagination
                page={page}
                totalCount={
                  enrichmentResults?.background?.enrich?.totalCount
                    ? enrichmentResults?.background?.enrich.totalCount
                    : undefined
                }
                pageSize={pageSize}
                onChange={(page) => {
                  setQueryString({ page: `${page}`, q: term, sortBy: sortBy, sortByDir: sortByDir });
                }}
              />
            </div>
          ) : null}
        </>
      ) : (
        <></>
      )}
      {tab === 2 && enrichmentResults ? (
        <TermVis
          enrichedTerms={
            termEnrichmentResults?.enrichedFunctionalTerms as EnrichedTermResult[]
          }
          sourceType={sourceType}
          setSourceType={setSourceType}
          setTab={setTab}
          setFilterTerm={setQueryString}
        />
      ) : (
        <></>
      )}
      {tab === 3 && enrichmentResults ? (
        <EnrichrTermVis
          enrichedTerms={
            EnrichrTermEnrichmentResults?.enrichedEnrichrTerms as EnrichrResult[]
          }
          setFilterTerm={setQueryString}
          setTab={setTab}
        />
      ) : (
        <></>
      )}
    </div>
  );
}

function GeneSetModalWrapper(props: {
  modalGeneSet: GeneSetModalT;
  setModalGeneSet: React.Dispatch<React.SetStateAction<GeneSetModalT>>;
}) {
  const { data: geneSet } = useViewGeneSetQuery({
    skip: props.modalGeneSet?.type !== "GeneSet",
    variables:
      props.modalGeneSet?.type === "GeneSet"
        ? {
            id: props.modalGeneSet.id,
          }
        : undefined,
  });
  const { data: overlap } = useOverlapQueryQuery({
    skip: props.modalGeneSet?.type !== "GeneSetOverlap",
    variables:
      props.modalGeneSet?.type === "GeneSetOverlap"
        ? {
            id: props.modalGeneSet.id,
            genes: props.modalGeneSet?.genes,
          }
        : undefined,
  });
  return (
    <GeneSetModal
      showModal={props.modalGeneSet !== undefined}
      term={props.modalGeneSet?.description}
      geneset={
        props.modalGeneSet?.type === "GeneSet"
          ? geneSet?.geneSet?.genes.nodes.map((gene) => gene.symbol)
          : props.modalGeneSet?.type === "GeneSetOverlap"
          ? overlap?.geneSet?.overlap.nodes.map((gene) => gene.symbol)
          : props.modalGeneSet?.type === "UserGeneSet"
          ? props.modalGeneSet.genes
          : undefined
      }
      setShowModal={(show) => {
        if (!show) props.setModalGeneSet(undefined);
      }}
    />
  );
}

function EnrichrModalWrapper(props: {
  modalGeneSet: GeneSetModalT;
  setModalGeneSet: React.Dispatch<React.SetStateAction<GeneSetModalT>>;
}) {
  const { data: geneSet } = useViewGeneSetQuery({
    skip: props.modalGeneSet?.type !== "GeneSet",
    variables:
      props.modalGeneSet?.type === "GeneSet"
        ? {
            id: props.modalGeneSet.id,
          }
        : undefined,
  });
  return (
    <EnrichrModal
      term={props.modalGeneSet?.description}
      geneset={geneSet?.geneSet?.genes.nodes.map((gene) => gene.symbol)}
      showModal={props.modalGeneSet !== undefined}
      setShowModal={(show) => {
        if (!show) props.setModalGeneSet(undefined);
      }}
    />
  );
}

function HypothesisModalWrapper(props: {
  modalGeneSet: GeneSetModalT;
  setModalGeneSet: React.Dispatch<React.SetStateAction<GeneSetModalT>>;
}) {
  const { data: overlap } = useOverlapQueryQuery({
    skip: props.modalGeneSet?.type !== "GeneSetOverlap",
    variables:
      props.modalGeneSet?.type === "GeneSetOverlap"
        ? {
            id: props.modalGeneSet.id,
            genes: props.modalGeneSet?.genes,
          }
        : undefined,
  });
  return (
    <HypothesisModal
      showModal={props.modalGeneSet !== undefined}
      term={props.modalGeneSet?.description}
      geneset={overlap?.geneSet?.overlap.nodes.map((gene) => gene.symbol)}
      gseId={
        props.modalGeneSet?.type === "GeneSetOverlap"
          ? props.modalGeneSet.gseId
          : undefined
      }
      setShowModal={(show) => {
        if (!show) props.setModalGeneSet(undefined);
      }}
    />
  );
}

function SamplesModalWrapper(props: {
  samples: string[];
  condition: string;
  setModalSamples: React.Dispatch<React.SetStateAction<string[] | undefined>>;
}) {
  return (
    <SamplesModal
      samples={props.samples}
      showModal={props.samples.length != 0 && props.condition !== ""}
      condition={props.condition}
      setShowModal={(show) => {
        if (!show) props.setModalSamples([]);
      }}
    />
  );
}

export default function Enrich({
  searchParams,
}: {
  searchParams: {
    dataset: string | string[] | undefined;
  };
}) {
  const dataset = ensureArray(searchParams.dataset)[0];
  const { data: userGeneSet } = useFetchUserGeneSetQuery({
    skip: !dataset,
    variables: { id: dataset },
  });
  const [modalGeneSet, setModalGeneSet] = React.useState<GeneSetModalT>();
  const [modalHypothesis, setModalHypothesis] = React.useState<GeneSetModalT>();
  const [modalEnrichr, setModalEnrichr] = React.useState<GeneSetModalT>();
  const [modalSamples, setModalSamples] = React.useState<string[]>();
  const [modalCondition, setModalCondition] = React.useState<string>();

  return (
    <>
      <div className="flex flex-row gap-2 alert bg-neutral-900 bg-opacity-10">
        <span className="prose">Input:</span>
        <label
          htmlFor="geneSetModal"
          className="prose underline cursor-pointer"
          onClick={(evt) => {
            setModalGeneSet({
              type: "UserGeneSet",
              genes: (userGeneSet?.userGeneSet?.genes ?? []).filter(
                (gene): gene is string => !!gene
              ),
              description: userGeneSet?.userGeneSet?.description || "Gene set",
            });
          }}
        >
          {userGeneSet?.userGeneSet?.description || "Gene set"}
          {userGeneSet ? (
            <> ({userGeneSet?.userGeneSet?.genes?.length ?? "?"} genes)</>
          ) : null}
        </label>
      </div>
      <EnrichmentResults
        userGeneSet={userGeneSet}
        setModalGeneSet={setModalGeneSet}
        setModalSamples={setModalSamples}
        setModalCondition={setModalCondition}
        setModalHypothesis={setModalHypothesis}
        setModalEnrichr={setModalEnrichr}
      />
      <SamplesModalWrapper
        samples={modalSamples ?? []}
        condition={modalCondition ?? ""}
        setModalSamples={setModalSamples}
      />
      <GeneSetModalWrapper
        modalGeneSet={modalGeneSet}
        setModalGeneSet={setModalGeneSet}
      />
      <HypothesisModalWrapper
        modalGeneSet={modalHypothesis}
        setModalGeneSet={setModalHypothesis}
      />
      <EnrichrModalWrapper
        modalGeneSet={modalEnrichr}
        setModalGeneSet={setModalEnrichr}
      />
    </>
  );
}
