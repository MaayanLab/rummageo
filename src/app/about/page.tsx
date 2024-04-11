import Image from "next/image";
import Stats from "../stats";
import Link from "next/link";

export default function About() {
  return (
    <div className="prose">
      <h2 className="title text-xl font-medium mb-3">Abstract</h2>
      <div className="flex">
        <div className="flex-col">
          <Image
            className={"rounded float-right ml-5"}
            src={"/images/rummageo_logo.png"}
            width={250}
            height={250}
            alt={"Rummageo"}
          ></Image>
          <p className="text-justify">
            The Gene Expression Omnibus (GEO) is a major open biomedical
            research repository for transcriptomics and other omics datasets. It
            currently contains millions of gene expression samples from tens of
            thousands of studies collected by many biomedical research
            laboratories from around the world. While users of the GEO
            repository can search the metadata describing studies and samples
            for locating relevant studies, there is currently no method or
            resource that facilitates global search of GEO at the data level. To
            address this shortcoming, we developed RummaGEO, a webserver
            application that enables gene expression signature search against
            all human and mouse RNA-seq studies deposited into GEO. To enable
            such a search engine, we performed offline automatic identification
            of conditions from uniformly aligned GEO studies available from
            ARCHS4, and then computed differential expression signatures to
            extract gene sets from these signatures. In total, RummaGEO
            currently contains <Stats show_human_gene_sets bold /> and{" "}
            <Stats show_mouse_gene_sets bold /> from <Stats show_gses bold />.
            Overall, RummaGEO provides an unprecedented resource for the
            biomedical research community enabling hypotheses generation for
            many future studies.
          </p>
          <h2 className="title text-xl font-medium mb-3 mt-3">Methods</h2>
          <p className="text-justify">
            We considered any GEO study aligned by{" "}
            <Link
              href="https://maayanlab.cloud/archs4/"
              className="underline cursor-pointer"
              target="_blank"
            >
              ARCHS4
            </Link>{" "}
            with at least three samples per condition with at least six samples
            in total collected for the study. Studies with more than 50 samples
            were discarded because such studied typically contain patient data
            that is not amenable for simple signature computation that compares
            two conditions. Samples were grouped using metadata provided by the
            GEO study. Specifically, K-means clustering of the embedding of
            concatenated sample{" "}
            <span style={{ fontStyle: "italic" }}>title</span>,{" "}
            <span style={{ fontStyle: "italic" }}>characteristic_ch1</span>, and{" "}
            <span style={{ fontStyle: "italic" }}>source_ch1 </span>
            fields were used to classify conditions. To create condition titles,
            common words across all samples for each condition were retained.
            Limma voom was used to compute differential expression signatures
            for each condition against all other conditions within each study.
            Additionally, we attempted to first identify any control conditions
            based on metadata and a discrete list of keywords that describe
            control conditions, for example, “wildtype”, “ctrl”, or “DSMO”. If
            such terms were identified, they were used to compare to the samples
            labeled with such term to all other condition groups. Up and down
            gene sets were extracted from each signature for genes with an
            adjusted p-value of less than 0.05. If less than five genes met this
            threshold, the gene set was discarded. If more than 2000 genes met
            this threshold, the threshold was lowered incrementally to 0.05,
            0.01, 0.005, and 0.001, until less than 2000 genes were retained.
            Additionally, we calculate a data-level confidence score of the
            condition groups with a silhouette score based on a PCA of
            normalized expression data wherein a value of 1 indicates perfect
            clustering and -1 indicates poor clustering.
          </p>
          <br></br>
          <p>
            For for information about using RummaGEO, please refer to the {" "}
            <Link
              href="/usermanual"
              className="underline cursor-pointer"
            >
              User Manual
            </Link>
            .
          </p>
          <br />
          <p>
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
          <p>
            This site is programatically accessible via a{" "}
            <Link
              href="/graphiql"
              className="underline cursor-pointer"
              target="_blank"
            >
              GraphQL API
            </Link>
            .
          </p>
          <br />
          <p>
            ReummaGEO is actively being developed by{" "}
            <Link
              className="underline cursor"
              href="https://labs.icahn.mssm.edu/maayanlab/"
              target="_blank"
            >
              the Ma&apos;ayan Lab
            </Link>
          </p>
        </div>
      </div>
    </div>
  );
}
