import React from "react";
import classNames from "classnames";
import Loading from "./loading";
import { useQueryGseSummaryQuery } from "@/graphql";
import parse from 'html-react-parser';
import fetchHypothesis from "@/utils/fetchHypothesis";
import getEnrichrTerms from "@/utils/getEnrichrTerms";

export default function HypothesisModal({
  geneset,
  term,
  gseId,
  showModal,
  setShowModal,
}: {
  geneset?: (string | null)[] | undefined;
  term: string | null | undefined;
  showModal?: boolean;
  gseId?: string | undefined;
  setShowModal: (show: boolean) => void;
}) {
  const [hypothesis, setHypothesis] = React.useState<Record<string, string> | null>(null);
  const [enrichrTerms, setEnrichrTerms] = React.useState<Record<string, string[]> | null>(null);
  const [error, setError] = React.useState<string | null>(null);
  const [loading, setLoading] = React.useState<boolean>(false);
  const [geneSetDesc, setGeneSetDesc] = React.useState<string>("");

  const wordCount = React.useMemo(
    () => geneSetDesc.split(/\s+/).filter((word) => word.length > 0).length,
    [geneSetDesc]
  );

  const { data } = useQueryGseSummaryQuery({
    variables: {
      id: gseId,
    },
  });

  const generateHypothesis = React.useCallback(async (t: string) => {
    if (wordCount < 100) {
      setError("The gene set description must be at least 100 words");
      setTimeout(() => setError(null), 3000);
      return;
    } else if (wordCount > 250) {
      setError("The gene set description can be at most 250 words");
      setTimeout(() => setError(null), 3000);
      return;
    } else if (t === '') {
        return
    }
    setLoading(true);
    try {
        const enrich = await getEnrichrTerms(geneset as string[], ['WikiPathway_2023_Human', 'GWAS_Catalog_2023', 'GO_Biological_Process_2023', 'MGI_Mammalian_Phenotype_Level_4_2021'])
        if (!enrich) {
            setLoading(false);
            setError('Error fetching enrichment results');
            return;
          } else {
        console.log(enrich)
        const hypothesisRes = await fetchHypothesis(
            geneSetDesc,
            data?.gseInfo?.summary || "",
            t || "",
            enrich[0] as Record<string, string[]>,
            enrich[1] as Record<string, string[]>
          );
          console.log(hypothesisRes)
          const newHypothesis = {...hypothesis}
          if (hypothesisRes != undefined) newHypothesis[t || ''] = hypothesisRes
          console.log(Object.keys(newHypothesis))
          setHypothesis(newHypothesis);
          setLoading(false);
        }
    }
    catch (e) {
      console.error(e);
      setError(e + '. Please try again.');
      setLoading(false);
      return;
    }
  }, [geneSetDesc, wordCount, data?.gseInfo?.summary, geneset, hypothesis]);

  return (
    <div
      className="z-40"
      onClick={(e) => {
        const hypModal = document.getElementById("hypothesis-div");
        if (!hypModal) return;
        const rect = hypModal?.getBoundingClientRect();

        const clickedInDialog =
          rect.top <= e.clientY &&
          e.clientY <= rect.top + rect.height &&
          rect.left <= e.clientX &&
          e.clientX <= rect.left + rect.width;
        console.log(clickedInDialog);
        if (!clickedInDialog) setShowModal(false);
      }}
    >
      {showModal ? (
        <>
          <div className="justify-center items-center flex overflow-x-hidden overflow-y-scroll fixed inset-0 z-30 focus:outline-none">
            <div className="relative w-auto my-6 mx-auto max-w-3xl">
              <div
                id="hypothesis-div"
                className="border-0 rounded-lg shadow-lg relative flex flex-col w-full bg-white outline-none focus:outline-none dark:bg-neutral-900"
              >
                <div className="p-3 border-b">
                  <p className="text-md text-center text-gray-900 dark:text-white">
                    Gene Set Overlap (
                    {geneset ? geneset?.length : "n"})
                  </p>
                </div>
                <p className="m-4 p-2">
                  {(hypothesis && Object.keys(hypothesis || {}).includes(term || '')) ? <>Hypothesis for <i>{term?.replace(".tsv", "")}</i> </>:<>Generate a hypothesis using GPT-4 to help speculate why a
                  signifigant overlap was found between your entered gene set
                  and the signature: <i>{term?.replace(".tsv", "")}</i>
                  <br></br>
                  <p className="mt-2">We will use the abstract of the GEO gene set and your submitted description along with highly enriched terms from Enrichr to generate a hypothesis.</p>
                  </>}
                  
                </p>
                {(hypothesis && Object.keys(hypothesis || {}).includes(term || '')) ? (
                  <div className="flex flex-col text-center justify-center mx-auto">
                    <div className="p-5 m-5 mt-0 text-left border-2 border-slate-300 rounded-lg font-light max-h-96 overflow-y-scroll break-all">{parse(hypothesis[term || ''])}</div>
                    <button
                      className="btn btn-sm btn-outline text-xs p-2 m-2"
                      type="button"
                      onClick={() => {
                        navigator.clipboard.writeText(hypothesis[term || ''].replace(/<[^>]+>/g, ''));
                      }}
                    >
                      Copy to Clipboard
                    </button>
                    <p className="font-extralight text-sm m-2">*Please use caution when interpreting LLM generated hypotheses*</p>
                  </div>
                ) : (
                  <></>
                )}

                {loading ? (
                  <Loading></Loading>
                ) : (<> {(hypothesis && Object.keys(hypothesis || {}).includes(term || '')) ? null :
                  <div className="flex flex-col">
                    <p className="m-4 mb-1 p-2">
                      To begin please enter a description of your gene set
                      between 100 and 250 words:
                    </p>

                    <textarea
                      value={geneSetDesc}
                      onChange={(evt) => {
                        setGeneSetDesc(evt.currentTarget.value);
                      }}
                      rows={8}
                      className="justify-center mx-auto textarea textarea-bordered w-3/4"
                      placeholder="Enter a description of your gene set here..."
                    />
                    <p className="m-0 p-1 text-center text-sm font-extralight">
                      {wordCount} words
                    </p>

                    <div className="flex items-center justify-center p-6 border-t border-solid border-slate-200 rounded-b">
                      <button
                        className="btn btn-sm btn-outline text-xs p-2 m-2"
                        onClick={() => generateHypothesis(term || '')}
                      >
                        Generate Hypothesis
                      </button>
                     
                    </div>
                  </div>}
                  </>
                )}
                <div className="flex flex-col justify-center text-center mx-auto">
                  <span
                    className={classNames("loading", "w-6", {
                      hidden: !loading,
                    })}
                  ></span>
                  <div
                    className={classNames("alert alert-error w-fit m-2 mt-0", {
                      hidden: !error,
                    })}
                  >
                    {error ?? null}
                  </div>
                </div>
              </div>
            </div>
          </div>
        </>
      ) : null}
    </div>
  );
}