import { EnrichedTermResult } from "@/graphql";
import { useState, useEffect } from "react";
import Pagination from "@/components/pagination";
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend,
} from "chart.js";
import { Bar } from "react-chartjs-2";
import WordCloud from "react-d3-cloud";
ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  Title,
  Tooltip,
  Legend
);
ChartJS.defaults.borderColor = '#85ffd6';
ChartJS.defaults.color = '#ffffff';


export const options = {
  responsive: true,
  indexAxis: "y" as const,
  elements: {
    bar: {
      borderWidth: 2,
    },
  },
  plugins: {
    legend: {
      position: "top" as const,
      labels: {
        color: "white",
      }
    },
    title: {
      display: true,
      text: "Top 10 Enriched Terms",
    },
  },
};

export default function TermVis({
  enrichedTerms,
}: {
  enrichedTerms: EnrichedTermResult[] | null | undefined;
}) {
  const [currentPage, setCurrentPage] = useState(1);
  const entriesPerPage = 10;
  const startIndex = (currentPage - 1) * entriesPerPage;
  const endIndex = startIndex + entriesPerPage;

  const totalPages = Math.ceil(enrichedTerms?.length || 1 / entriesPerPage);

  const barChartData = {
    labels: enrichedTerms?.slice(startIndex, endIndex).map((t) => t.term) || [],
    datasets: [
      {
        label: "-log10(adjPvalue)",
        data:
          enrichedTerms
            ?.slice(startIndex, endIndex)
            .map((t) => -Math.log10(t?.adjPvalue || 0)) || [],
        backgroundColor: "rgba(75, 192, 192, 0.2)",
        borderColor: "rgba(75, 192, 192, 1)",
        borderWidth: 2,
      },
    ],
  };

  // Prepare data for word cloud
  const wordCloudData =
    enrichedTerms?.map((term) => ({
      text: term.term || "",
      value: term.count || 0,
    })) || [];

  return (
    <div className="flex flex-col w-full">
      <div className="mx-auto min-h-full h-fit w-10/12">
      <WordCloud data={wordCloudData} rotate={0} padding={0} width={300} height={200}/>
            <Bar data={barChartData} options={options} />
            
      </div>
      <table className="table table-xs">
        <thead>
          <tr>
            <th>Term</th>
            <th>Count</th>
            <th>Odds Ratio</th>
            <th>P-value</th>
            <th>Adj. P-value</th>
          </tr>
        </thead>
        <tbody>
          {enrichedTerms?.slice(startIndex, endIndex).flatMap((row, i) => {
            return (
              <tr key={i}>
                <td>
                  <b>{row.term}</b>
                </td>
                <td>{row.count}</td>
                <td>{row.oddsRatio?.toPrecision(3)}</td>
                <td>{row.pvalue?.toPrecision(3)}</td>
                <td>{row.adjPvalue?.toPrecision(3)}</td>
              </tr>
            );
          })}
        </tbody>
      </table>
      <div className="w-full flex flex-col items-center mt-2">
        <Pagination
          page={currentPage}
          totalCount={totalPages}
          pageSize={entriesPerPage}
          onChange={setCurrentPage}
        />
      </div>
    </div>
  );
}
