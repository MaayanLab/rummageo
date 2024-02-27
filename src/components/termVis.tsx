import { EnrichedTermResult } from "@/graphql";
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

export default function TermVis({
  enrichedTerms,
}: {
  enrichedTerms: EnrichedTermResult[] | null | undefined;
}) {
  const [currentPage, setCurrentPage] = useState(1);
  const entriesPerPage = 10;
  const startIndex = (currentPage - 1) * entriesPerPage;
  const endIndex = startIndex + entriesPerPage;
  const [searchTerm, setSearchTerm] = useState("");

  

  const enrichedTermsFiltered = useMemo(() => enrichedTerms?.filter(r => r.term?.toLowerCase().includes(searchTerm.toLowerCase())), [enrichedTerms, searchTerm]);

  const totalPages = useMemo(() => Math.ceil(enrichedTermsFiltered?.length || 1 / entriesPerPage), [enrichedTermsFiltered, entriesPerPage])

  const barChartData = {
    labels: enrichedTermsFiltered?.slice(startIndex, endIndex).map((t) => t.term) || [],
    datasets: [
      {
        label: "-log10(adjPvalue)",
        data:
          enrichedTermsFiltered
            ?.slice(startIndex, endIndex)
            .map((t) => -Math.log10(t?.adjPvalue || 0)) || [],
        backgroundColor: "#13ad65",
        borderColor: "#fcfffe",
        borderWidth: 2,
      },
    ],
  };

  const wordCloudData = useMemo(() => enrichedTerms?.map((term) => ({
    text: term.term || "",
    value: (term.count ?? 1) * 10,
  })) || [], [enrichedTerms]);

  const chartRef = useRef<ChartJS<"bar", (number | [number, number] | Point | BubbleDataPoint | null)[], unknown> | null>(null);

  return (
    <div className="flex flex-col w-full">
      <div className="mx-auto min-h-full h-fit w-10/12">
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
                  const blob = blobTsv(['term', 'CountInEnrichedGSEs', 'OtherTermsInEnrichedGSEs', 'TotalGSEsWithTerm', 'TotalOtherTerms', 'OddsRatio', 'Pvalue', 'AdjPvalue'], enrichedTermsFiltered, elt => {
                  return {
                  term: elt.term,
                  CountInEnrichedGSEs: elt.count,
                  OtherTermsInEnrichedGSEs: elt.notTermCount,
                  TotalGSEsWithTerm: elt.totalTermCount,
                  TotalOtherTerms: elt.totalNotTermCount,
                  OddsRatio: elt.oddsRatio,
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
            <th>Count in Enriched GSEs</th>
            <th>Other terms in Enriched GSEs</th>
            <th>Total GSEs w/ term</th>
            <th>Total other terms</th>
            <th>Odds Ratio</th>
            <th>P-value</th>
            <th>Adj. P-value</th>
          </tr>
        </thead>
        <tbody>
          {enrichedTermsFiltered?.slice(startIndex, endIndex).flatMap((row, i) => {
            return (
              <tr key={i} className="text-center">
                <td>
                  <b>{row.term}</b>
                </td>
                <td>{row.count}</td>
                <td>{row.notTermCount}</td>
                <td>{row.totalTermCount}</td>
                <td>{row.totalNotTermCount}</td>
                <td>{row.oddsRatio?.toPrecision(3)}</td>
                <td>{row.pvalue?.toPrecision(3)}</td>
                <td>{row.adjPvalue?.toPrecision(3)}</td>
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
        <div className="inline-flex w-fit">
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
  );
}
