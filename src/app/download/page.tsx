'use client'
import React from 'react'
import Stats from "@/app/stats";
import { useLatestReleaseQuery } from "@/graphql";

export default function Download() {
  const { data } = useLatestReleaseQuery()
  const latest_release_date = React.useMemo(() => {
    const date = new Date(data?.releases?.nodes[0]?.created)
    return date.toDateString()
  }, [data])
  return (
    <div className="prose">
      <h2 className="title text-xl font-medium mb-3">Downloads</h2>
      <br />
      <p>
        Gene set libraries are available for download in GMT format. Accompanying metadata is available in JSON format which contains gene set condition titles and GSM ids.
      </p>
      <div className="grid lg:grid-cols-2 gap-4 my-4">
        <a className="stats shadow" href="https://s3.amazonaws.com/maayanlab-public/rummageo/human_12-04-2023.gmt.gz" download="human-geo-auto.gmt.gz">
          <div className="stat gap-2">
            <div className="stat-title">human-geo-auto.gmt.gz</div>
            <div className="stat-value text-xl">135,264 gene sets</div>
            <div className="stat-desc whitespace-normal">
            <span className="whitespace-nowrap">386MB compressed</span>, <span className="whitespace-nowrap">874MB uncompressed</span>, <span className="whitespace-nowrap">Last Updated Fri Dec 8 2023</span>
            </div>
          </div>
        </a>
        <a className="stats shadow" href="https://s3.amazonaws.com/maayanlab-public/rummageo/gse_processed_meta_human.json" download="geo_processed_meta_human.json">
          <div className="stat gap-2">
            <div className="stat-title">geo_processed_meta_human.json</div>
            <div className="stat-value text-xl">human metadata</div>
            <div className="stat-desc whitespace-normal">
              <span className="whitespace-nowrap">11.7MB uncompressed</span>, <span className="whitespace-nowrap">Last Updated Fri Dec 8 2023</span>
            </div>
          </div>
        </a>
        <a className="stats shadow" href="https://s3.amazonaws.com/maayanlab-public/rummageo/mouse_12-04-2023.gmt.gz" download="mouse-geo-auto.gmt.gz">
          <div className="stat gap-2">
            <div className="stat-title">mouse-geo-auto.gmt.gz</div>
            <div className="stat-value text-xl">158,062 gene sets</div>
            <div className="stat-desc whitespace-normal">
              <span className="whitespace-nowrap">418MB compressed</span>, <span className="whitespace-nowrap">878MB uncompressed</span>, <span className="whitespace-nowrap">Last Updated Fri Dec 8 2023</span>
            </div>
          </div>
        </a>
        <a className="stats shadow" href="https://s3.amazonaws.com/maayanlab-public/rummageo/gse_processed_meta_mouse.json" download="gse_processed_meta_mouse.json">
          <div className="stat gap-2">
            <div className="stat-title">gse_processed_meta_mouse.json</div>
            <div className="stat-value text-xl">mouse metadata</div>
            <div className="stat-desc whitespace-normal">
              <span className="whitespace-nowrap">13.6MB</span>, <span className="whitespace-nowrap">Last Updated Fri Dec 8 2023</span>
            </div>
          </div>
        </a>
        <a className="stats shadow" href="https://s3.amazonaws.com/maayanlab-public/rummageo/random50k_co-occurrence_coding.f" download="human_random50k_co-occurrence_coding.f">
          <div className="stat gap-2">
            <div className="stat-title">human_random50k_co-occurrence_coding.f</div>
            <div className="stat-value text-xl">human coding co-occurence</div>
            <div className="stat-desc whitespace-normal">
              <span className="whitespace-nowrap">1.13G</span>, <span className="whitespace-nowrap">Last Updated Thu Apr 4 2024</span>
            </div>
          </div>
        </a>
        <a className="stats shadow" href="https://s3.amazonaws.com/maayanlab-public/rummageo/random50k_co-occurrence_noncoding.f" download="human_random50k_co-occurrence_noncoding.f">
          <div className="stat gap-2">
            <div className="stat-title">human_random50k_co-occurrence_noncoding.f</div>
            <div className="stat-value text-xl">human noncoding co-occurence</div>
            <div className="stat-desc whitespace-normal">
              <span className="whitespace-nowrap">2.41G</span>, <span className="whitespace-nowrap">Last Updated Thu Apr 4 2024</span>
            </div>
          </div>
        </a>
        <a className="stats shadow" href="https://s3.amazonaws.com/maayanlab-public/rummageo/mouse_random50k_co-occurrence_coding.f" download="mouse_random50k_co-occurrence_coding.f">
          <div className="stat gap-2">
            <div className="stat-title">mouse_random50k_co-occurrence_coding.f</div>
            <div className="stat-value text-xl">mouse coding co-occurence</div>
            <div className="stat-desc whitespace-normal">
              <span className="whitespace-nowrap">1.16G</span>, <span className="whitespace-nowrap">Last Updated Thu Apr 4 2024</span>
            </div>
          </div>
        </a>
        <a className="stats shadow" href="https://s3.amazonaws.com/maayanlab-public/rummageo/mouse_random50k_co-occurrence_noncoding.f" download="mouse_random50k_co-occurrence_noncoding.f">
          <div className="stat gap-2">
            <div className="stat-title">mouse_random50k_co-occurrence_noncoding.f</div>
            <div className="stat-value  text-xl">mouse noncoding co-occurence</div>
            <div className="stat-desc whitespace-normal">
              <span className="whitespace-nowrap">1.41G</span>, <span className="whitespace-nowrap">Last Updated Thu Apr 4 2024</span>
            </div>
          </div>
        </a>
      </div>
      <p>
        Developed in <a className='underline cursor' href="https://labs.icahn.mssm.edu/maayanlab/">the Ma&apos;ayan Lab</a>
      </p>
    </div>
  )
}
