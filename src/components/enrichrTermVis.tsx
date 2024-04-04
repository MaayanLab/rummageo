import { EnrichrResult } from "@/graphql";
import { useState, useRef, useMemo } from "react";
import Pagination from "@/components/pagination";
import clientDownloadBlob from '@/utils/clientDownloadBlob';
import blobTsv from '@/utils/blobTsv';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend,
  Point,
  BubbleDataPoint,
} from "chart.js";
import { Bar } from "react-chartjs-2";
import { FaDownload } from "react-icons/fa";
import { Wordcloud } from "@visx/wordcloud";
import { Text } from "@visx/text";
import { FaSearch } from "react-icons/fa";
import { TiDeleteOutline } from "react-icons/ti";
import { MdOutlineFileDownload } from "react-icons/md";
import classNames from "classnames";
import Loading from "@/components/loading";
const d3ToPng = require('d3-svg-to-png');


ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
);
ChartJS.defaults.borderColor = "#ffffff";
ChartJS.defaults.color = "#518272";
ChartJS.defaults.backgroundColor = "#000000";

const colors = ["#143059", "#2F6B9A", "#3c997a"];

export const options = {
  responsive: true,
  indexAxis: "y" as const,
  elements: {
    bar: {
      borderWidth: 1,
      color: "#00C453",
    },
  },
  plugins: {
    legend: {
      position: "bottom" as const,
      labels: {
        color: "#518272",
      },
    },
  },
};


export interface WordData {
  text: string;
  value: number;
}

export default function EnrichrTermVis({
  enrichedTerms,
  setFilterTerm,
  setTab
}: {
  enrichedTerms: EnrichrResult[] | undefined,
  setFilterTerm: Function,
  setTab: React.Dispatch<React.SetStateAction<number>>,
}) {
  const [currentPage, setCurrentPage] = useState(1);
  const entriesPerPage = 10;
  const startIndex = (currentPage - 1) * entriesPerPage;
  const endIndex = startIndex + entriesPerPage;
  const [searchTerm, setSearchTerm] = useState("");

  const enrichedTermsFiltered = useMemo(() => enrichedTerms?.filter(r => r?.term?.toLowerCase().includes(searchTerm.toLowerCase())), [enrichedTerms, searchTerm]);

  const totalPages = useMemo(() => Math.ceil(enrichedTermsFiltered?.length || 1 / entriesPerPage), [enrichedTermsFiltered, entriesPerPage])

  const barChartData = {
    labels: enrichedTermsFiltered?.slice(startIndex, endIndex).map((t) => t?.term) || [],
    datasets: [
      {
        label: "-log10(P-value)",
        data:
          enrichedTermsFiltered
            ?.slice(startIndex, endIndex)
            .map((t) => -Math.log10(t?.pvalue || 0)) || [],
        backgroundColor: "#13ad65",
        borderColor: "#fcfffe",
        borderWidth: 2,
      },
    ],
  };

  const wordCloudData = useMemo(() => enrichedTerms?.map((term) => ({
    text: term.term || "",
    value: Math.max((-Math.log(term?.pvalue ?? .5) * 50), 200),
  })) || [], [enrichedTerms]);

  const chartRef = useRef<ChartJS<"bar", (number | [number, number] | Point | BubbleDataPoint | null)[], unknown> | null>(null);

  const handleDownloadImage = async () => {
    const elt = document.getElementById('wc-svg')?.children[0] as HTMLElement
    d3ToPng(elt, 'enriched_terms_wc', {download: true, format: 'png', ignore: '.ignore',  quality: 1})
  };


  return (
    <>
    {enrichedTerms ? 
    <div className="flex flex-col w-full">
      <div className="mx-auto min-h-full h-fit w-10/12">
        <p className="font-light mb-3">Top 3 significant (adj. p-value &lt; 0.05) <a href="https://maayanlab.cloud/Enrichr/" target="_blank">Enrichr</a> terms precomputed for 
        each signature appearing in top 500 signatures from selected libraries (ChEA_2022, KEGG_2021_Human, 
        WikiPathway_2023_Human, MGI_Mammalian_Phenotype_Level_4_2021, GO_Biological_Process_2023, Human_Phenotype_Ontology, GWAS_Catalog_2023). Signfigance of term rankings is computed with the Kolmogorov–Smirnov test.</p>
        <div className="flex gap-3">
        </div>
        <button
          className="float-right m-3"
          onClick={(evt) => {
            evt.preventDefault();
            if (chartRef.current != null) {
              var a = document.createElement("a");
              a.href = chartRef.current.toBase64Image();
              a.download = "enriched_terms.png";
              a.click();
            }
          }}
        >
          <FaDownload />
        </button>

        <Bar data={barChartData} options={options} ref={chartRef} />
      </div>
      <form
          className="join flex flex-row place-content-end place-items-center"
          onSubmit={(evt) => {
            evt.preventDefault();
          }}
        >
          <input
            type="text"
            className="input input-bordered bg-transparent join-item"
            value={searchTerm}
            onChange={(evt) => {
              setSearchTerm(evt.currentTarget.value);
              setCurrentPage(1);
            }}
          />
          <div className="tooltip" data-tip="Search results">
            <button type="submit" className="btn join-item bg-transparent ml-2">
              <FaSearch />
            </button>
          </div>
          <div className="tooltip" data-tip="Clear search">
            <button
              type="reset"
              className="btn join-item bg-transparent"
              onClick={(evt) => {
                evt.preventDefault();
                setSearchTerm("");
              }}
            >
              <TiDeleteOutline />
            </button>
          </div>
          <a
            href={``}
            download="results.tsv"
          >
            <div className="tooltip" data-tip="Download results">
              <button
                type="button"
                className="btn join-item font-bold text-2xl pb-1 bg-transparent"
                onClick={(evt) => {
                  evt.preventDefault()
                  if (!enrichedTermsFiltered) return
                  const blob = blobTsv(['term', 'CountInEnrichedGSEs', 'Pvalue', 'AdjPvalue'], enrichedTermsFiltered, elt => {
                  return {
                  term: elt.term,
                  CountInEnrichedGSEs: elt.count,
                  Pvalue: elt.pvalue,
                  AdjPvalue: elt.adjPvalue
                  }
                })
                clientDownloadBlob(blob, 'results.tsv')
              }}
              >
                <MdOutlineFileDownload />
              </button>
            </div>
          </a>
        </form>
      <div className="overflow-x-scroll">
      <table className="table table-s border-2 border-white dark: mt-10">
        <thead>
          <tr className="text-center">
            <th>Term</th>
            <th>P-value</th>
            <th>Adj. P-value</th>
            <th>Count in Enriched Signatures</th>
          </tr>
        </thead>
        <tbody>
          {enrichedTermsFiltered?.slice(startIndex, endIndex).flatMap((row, i) => {
            return (
              <tr key={i} className="text-center">
                <td>
                  <a className="link" onClick={() => {
                    setFilterTerm({page: "1", q: row.term})
                    setTab(1)
                  }
                  }>
                    <div className="tooltip tooltip-right underline" data-tip={`View enriched signatures which mention ${row.term}`}>
                      <b>{row.term}</b>
                    </div>
                  </a>
                </td>
                
                <td>{row.pvalue?.toPrecision(3)}</td>
                <td>{row.adjPvalue?.toPrecision(3)}</td>
                <td>{row.count}</td>
              </tr>
            );
          })}
        </tbody>
      </table>
      </div>
      <div className="w-full flex flex-col items-center mt-2 mb-10">
        <Pagination
          page={currentPage}
          totalCount={totalPages}
          pageSize={entriesPerPage}
          onChange={setCurrentPage}
        />
      </div>
      <div className="flex-col h-fit w-full mx-auto text-center justify-center">
      <button
          className="float-right m-3"
          onClick={handleDownloadImage}
        >
          <FaDownload/>
          </button>
        <div className="inline-flex w-fit">
        
        <div id="wc-svg">
          <Wordcloud
            
            words={wordCloudData}
            rotate={0}
            padding={5}
            width={800}
            height={200}
          >
            {(cloudWords) =>
              cloudWords.map((w, i) => (
                <Text
                  key={w.text}
                  fill={colors[i % colors.length]}
                  textAnchor={"middle"}
                  transform={`translate(${w.x}, ${w.y})`}
                  fontSize={w.size}
                  fontFamily={w.font}
                >
                  {w.text}
                </Text>
              ))
            }
          </Wordcloud>
          </div>
          
          </div>
      </div>
    </div> : <Loading />}
    </>
  );
}
