'use client';
import Image from "next/image";
import Link from "next/link";
import { FaCopy } from "react-icons/fa";

const codeExample = `import requests
import json

url = "https://rummageo.com/graphql"

def enrich_rummageo(geneset: list):
    query = {
    "operationName": "EnrichmentQuery",
        "variables": {
            "filterTerm": "",
            "offset": 0,
            "first": 100,
            "genes": geneset,
            "id": "15c56ba6-a293-4932-bcbc-27fc2e4327ab"
        },
        "query": """query EnrichmentQuery($genes: [String]!, $filterTerm: String = "", $offset: Int = 0, $first: Int = 10, $id: UUID!) {
            background(id: $id) {
                id
                species
                enrich(genes: $genes, filterTerm: $filterTerm, offset: $offset, first: $first) {
                nodes {
                    pvalue
                    adjPvalue
                    oddsRatio
                    nOverlap
                    geneSet {
                    id
                    term
                    nGeneIds
                    geneSetPmidsById {
                        nodes {
                        gse
                        gseId
                        pmid
                        sampleGroups
                        platform
                        publishedDate
                        title
                        __typename
                        }
                        __typename
                    }
                    __typename
                    }
                    __typename
                }
                totalCount
                __typename
                }
                __typename
            }
        }
        """
    }

    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json"
    }

    response = requests.post(url, data=json.dumps(query), headers=headers)

    if response.status_code == 200:
        res = response.json()
        return res
`;

export default function UserManual() {
  return (
    <div className="prose">
      <div className="flex">
        <div className="flex-col mx-auto max-w-5xl">
          <h2 className="title text-2xl font-medium mb-3 mt-3 text-center">
            RummaGEO Documentation
          </h2>
          <div className="navbar block text-center">
            <div className="navbar-center">
              <ul className="menu menu-horizontal gap-3 flex text-lg">
                <li>
                  <Link
                    href="/about"
                    className="underline cursor-pointer"
                    shallow
                  >
                    Abstract
                  </Link>
                </li>
                <li>
                  <Link
                    href="/about"
                    className="underline cursor-pointer"
                    shallow
                  >
                    Methods
                  </Link>
                </li>
                <li>
                  <Link
                    href="#gene-set-search"
                    className="underline cursor-pointer"
                    shallow
                  >
                    Gene Set Search
                  </Link>
                </li>
                <li>
                  <Link
                    href="#pubmed-search"
                    className="underline cursor-pointer"
                    shallow
                  >
                    PubMed Search
                  </Link>
                </li>
                <li>
                  <Link
                    href="#metadata-search"
                    className="underline cursor-pointer"
                    shallow
                  >
                    Metadata Search
                  </Link>
                </li>
                <li>
                  <Link
                    href="#api"
                    className="underline cursor-pointer"
                    shallow
                  >
                    API
                  </Link>
                </li>
              </ul>
            </div>
          </div>
          <br></br>
          <h2 className="title text-xl font-medium mb-3" id="gene-set-search">
            1.1 Gene Set Search
          </h2>
          <p>
          The Gene Set Search page enables users to search the RummaGEO database for gene sets that match their query gene set. Similarity to gene sets contained within the RummaGEO database with the query gene set is measured with Fisher&apos;s exact test. Any significantly overlapping gene sets are returned to the user along with their accompanying metadata. User query gene sets can be pasted or typed into the input form with each gene on a new line, or the user may upload a file containing genes where the genes are listed with new line, tab, or comma separators. Based on the gene symbols within the query gene set, this query will run against either the collection of automatically generated human gene sets or mouse gene sets:

          </p>
          <Image
            src="/images/gene-set-search-1.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <p>
          Paginated results are returned with the total number of gene sets which the query was compared to and the number of those gene sets which were significantly enriched. Enrichment statistics are provided on the right side of the table. The user may explore the metadata associated with the signatures on each results page, by inspecting the title of the GEO study, the corresponding linked GEO accession, and, if available, the linked PubMed ID:
          </p>
          <Image
            src="/images/gene-set-search-2.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <p>
          Additionally, the sample IDs (GSM) and their metadata associated with each condition are displayed in a modal box when clicked:
          </p>
          <Image
            src="/images/gene-set-search-3.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <p>
          Other information such as the directionality of the gene set, the platform (GPL), the date of the study, the number of overlapping genes, and the total number of genes in the enriched gene set are also displayed. Clicking the overlap or gene set size will open a modal box with the corresponding genes as well as buttons to copy the gene set to the clipboard, or to view the enrichment results on RummaGEO or Enrichr:
          </p>
          <Image
            src="/images/gene-set-search-4.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <p>
          To further filter and refine the results, the user may use the search bar located above the table to search for gene sets containing certain keywords. This allows, for instance, to view enriched results related to macrophages based upon the same input gene set. The total number of enriched gene sets is updated accordingly:
          </p>
          <Image
            src="/images/gene-set-search-5.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <p>
          Results may also be easily downloaded using the button to the far right of the search feature in a tab delimited format:
          </p>
          <Image
            src="/images/gene-set-search-6.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />

          <p>
          Selecting the Common Terms in Matching Gene Sets tab will provide enrichment results computed with the Fisher&apos;s Exact Test for the first 5,000 unique GSEs from signatures returned in the gene set search result:
          </p>
          <Image
            src="/images/gene-set-search-7.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
           <p>
          Additionally, below the bar chart, the user may view a table with all enriched terms including the signfigance of each term&apos;s appearance. Clicking a term will result in the user being brought back to the Matching Gene Sets tab wherein the user can view the gene sets that contain the selected term:
          </p>
          <Image
            src="/images/gene-set-search-8.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />

          <p>
         Selecting the Enrichr Terms tab will show the most commonly appearing Enrichr terms in the top 500 signatures returned from the gene set search. Enrichr terms are precomputed for all RummaGEO signatures for a selection of libraries (ChEA 2022, KEGG 2021 Human, WikiPathway 2023 Human, GO Biological Process 2023, MGI Mammalian Phenotype Level 4 2021, Human Phenotype Ontology, and GWAS Catalog 2023):
          </p>
          <Image
            src="/images/gene-set-search-9.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          
          
          <h3
            className="title text-lg font-medium mb-3 mt-10"
            id="pubmed-search"
          >
            Hypothesis Generation
          </h3>

          <p>
          Users may also generate hypotheses concerning their gene set&apos;s overlap with any RummaGEO signatures:
          </p>
          <Image
            src="/images/gene-set-search-10.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />

          <p>
          After clicking the generate hypothesis button, users will need to enter a description of their gene set, which will be provided to the LLM:
          </p>
          <Image
            src="/images/gene-set-search-11.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <p>
          RummaGEO takes this description together with the matching RummaGEO gene set study abstract, and the top three significantly enriched terms from the overlapping genes from several Enrichr libraries (WikiPathway 2023 Human, GWAS Catalog 2023, GO Biological Process 2023, MGI Mammalian Phenotype Level 4 2021). The prompt additionally instructs the large language model (LLM) to reference all the provided descriptions and context of the gene sets, as well as the highly enriched terms from Enrichr. Hypotheses are then parsed to find references to any enriched terms, and insert the enrichment statistics as part of the hypothesis description:

          </p>
          <Image
            src="/images/gene-set-search-12.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />


          <h2
            className="title text-xl font-medium mb-3 mt-10"
            id="pubmed-search"
          >
            1.2 PubMed Search
          </h2>
          <p>
          The PubMed Search page enables users to search for gene sets in RummaGEO based on a PubMed search using the PubMed API query. The top 5000 publications returned from the user&apos;s query are used to display extracted gene sets from the GEO studies associated with the returned papers. The number of articles returned by the PubMed API along with the number of associated gene sets and associated publications in the RummaGEO database are displayed at the top of the result:
          </p>
          <Image
            src="/images/pubmed-search-1.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <p>
          Paginated results are grouped by GEO study (GSE) with the corresponding signatures and available metadata located in a dropdown table:
          </p>
          <Image
            src="/images/pubmed-search-2.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <p>
          Additionally, the results can be further filtered using the search bar at the top right of the table and the as well as downloaded in a tab-delimited format:
          </p>
          <Image
            src="/images/pubmed-search-3.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <h2
            className="title text-xl font-medium mb-3 mt-10"
            id="metadata-search"
          >
            1.3 Metadata Search
          </h2>
          RummaGEO also provides direct metadata search of the GEO studies contained within the database. Paginated results are returned with accompanying metadata of the returned signatures:
          <Image
            src="/images/metadata-search-1.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
         These results can also be filtered using the search bar at the top right of the table and the results table can be downloaded in a tab-delimited format:

          <Image
            src="/images/metadata-search-2.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <h2 className="title text-xl font-medium mb-3 mt-10" id="api">
            1.4 API
          </h2>
          <p>
            {" "}
            RummaGEO priovides programtic access through a GraphQL endpoint.
            Users can learn more about GraphQL queries from their provided{" "}
            <a href="https://graphql.org/learn/" target="_blank">
              documentation
            </a>
            . The RummaGEO{" "}
            <Link href="/graphiql" target="_blank">
              GraphQL endpoint
            </Link>{" "}
            and asscoiated Postgres database provide users with a wide range of
            available queries and with a user interface to test and develop
            these queries:
          </p>
          <Image
            src="/images/api-1.png"
            width={600}
            height={500}
            alt={""}
            className="border rounded-lg mx-auto my-4"
          />
          <p>
          For example, enrichment analysis queries can be performed in Python agaisnt the RummaGEO human gene sets using the requests library as follows:
          </p>
          <div className="text-gray bg-slate-300 dark:bg-slate-700 text-xs font-mono mt-5 p-5 rounded-lg box-content sm:max-w-xl sm:overflow-scroll md:max-w-xl lg:max-w-3xl xl:max-w-full">
            <button className="float-right" onClick={() => navigator.clipboard.writeText(codeExample)}><FaCopy/></button>
            <pre>
              <code>
                {codeExample}
              </code>
            </pre>
          </div>
          <p className="mt-10">
            This database is updated with new releases of{" "}
            <Link
              href="https://maayanlab.cloud/archs4/index.html"
              className="underline cursor-pointer"
              target="_blank"
            >
              ARCHS4
            </Link>
            .
          </p>
          <br />
          <>
            RummaGEO is actively being developed by{" "}
            <Link
              className="underline cursor"
              href="https://labs.icahn.mssm.edu/maayanlab/"
              target="_blank"
            >
              the Ma&apos;ayan Lab
            </Link>
            .
          </>
        </div>
      </div>
    </div>
  );
}
