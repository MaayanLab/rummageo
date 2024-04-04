import React from "react";
import classNames from "classnames";
import Loading from "./loading";
import { useEnrichrStatsQuery } from "@/graphql";
import EnrichrButton from "./enrichrButton";

export default function EnrichrModal({
  geneset,
  term,
  showModal,
  setShowModal,
}: {
  geneset?: (string | null)[] | undefined;
  term: string | null | undefined;
  showModal?: boolean;
  setShowModal: (show: boolean) => void;
}) {
  var species = ''
  if (term?.includes('human')) {
    species = 'human'
  } else {
    species = 'mouse'
  }
  const [rows, setRows] = React.useState<any[]>([])
  const { data: data } = useEnrichrStatsQuery({
    variables: { sig: term || '', organism: species },
  })
  const enrichrStats = data?.enrichrTermBySigAndOrganism?.enrichrStats
  
  React.useEffect(() => {
    if (enrichrStats) {
      const rows: any[] = []
      Object.keys(enrichrStats).map((lib, i) => (
        enrichrStats[lib].map((row: any, j: any) => {
          if ( Number(row[2]) < 0.05) {
            rows.push({
              term: row[0],
              oddsRatio: Number(row[3]).toFixed(4),
              pValue: Number(row[1]).toExponential(2),
              adjPValue: Number(row[2]).toExponential(2),
              overlap: Number(row[4]),
              library: lib
            })
        }
      })
      ))
      rows.sort((a, b) => a.adjPValue - b.adjPValue);
      setRows(rows)
    }
  }, [enrichrStats])



  return (
    <div
      className="z-40"
      onClick={(e) => {
        const modal = document.getElementById("enrichr-div");
        if (!modal) return;
        const rect = modal?.getBoundingClientRect();

        const clickedInDialog =
          rect.top <= e.clientY &&
          e.clientY <= rect.top + rect.height &&
          rect.left <= e.clientX &&
          e.clientX <= rect.left + rect.width;
        if (!clickedInDialog) setShowModal(false);
      }}
    >
      {showModal ? (
        <>
          <div className="justify-center items-center flex overflow-x-hidden overflow-y-scroll fixed inset-0 z-30 focus:outline-none">
            <div className="relative w-auto my-6 mx-auto max-w-5xl ">
              <div
                id="enrichr-div"
                className="border-0 rounded-lg shadow-lg relative flex flex-col w-full bg-white outline-none focus:outline-none dark:bg-neutral-900 overflow-scroll "
              >
                <div className="p-3 border-b">
                  <p className="text-md text-center text-gray-900 dark:text-white">
                    Gene Set (
                    {geneset ? geneset?.length : "n"})
                  </p>
                </div>
                
                  {enrichrStats ? 
                  <>
                
                    <div className="flex flex-col justify-center text-center mx-auto ">
                      <table className="table table-xs">
                        <thead>
                          <tr>
                          <th>
                            Term
                          </th>
                          <th>
                            Odds Ratio
                          </th>
                          <th>
                            P-value
                          </th>
                          <th>
                            Adj. P-value
                          </th>
                          <th>
                            Overlap
                          </th>
                          <th>
                            Library
                          </th>
                          </tr>

                        </thead>
                        <tbody>
                          {rows.map((row, i) => 
                          
                          {
                            return (
                              <tr key={i}>
                                <td>
                                  <b>{row.term}</b>
                                </td>
                                <td>
                                  {Number(row.oddsRatio).toFixed(2)}
                                </td>
                                <td>
                                  {Number(row.pValue).toExponential(2)}
                                </td>
                                <td>
                                  {Number(row.adjPValue).toExponential(2)}
                                </td>
                                <td>
                                  {Number(row.overlap)}
                                </td>
                                <td>
                                  {row.library}
                                </td>
                              </tr>)}
                            )}
                        </tbody>
                      </table>
                    </div>
                  
                  </>
                   : <><Loading/></>}
                  <EnrichrButton genes={geneset} description={term}></EnrichrButton>
                <div className="flex flex-col justify-center text-center mx-auto">
                </div>
              </div>
          </div>
    </div></>) : <></>}
    </div>
  );
}