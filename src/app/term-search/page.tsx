'use client'
import React from 'react'
import { useSuspenseQuery } from '@apollo/experimental-nextjs-app-support/ssr'
import { GeneSetLibraryTermSearchDocument, GeneSetLibraryTermSearchQuery } from '@/graphql'
import ensureArray from '@/utils/ensureArray'
import { useRouter } from 'next/navigation'
import TermTable from '@/components/termTable'
import Loading from '@/components/loading'

function TermSearchResults({ terms }: { terms: string[] }) {
  const { data } = useSuspenseQuery<GeneSetLibraryTermSearchQuery>(GeneSetLibraryTermSearchDocument, {
    variables: {
      terms,
      first: 10000
    }
  })

  return (
    <ul>
      {data?.geneSetLibraries?.nodes
        .filter(geneSetLibrary => geneSetLibrary.termSearchCount.nodes.length > 0)
        .map((geneSetLibrary, i) => (
          <div key={geneSetLibrary.name} className='text-center mt-5'>
            <p className='text-lg'> Your search term is contained in {geneSetLibrary.termSearchCount.totalCount} gene sets.</p>
            <TermTable terms={geneSetLibrary.termSearchCount.nodes}></TermTable>
          </div>
        )) ?? null}
    </ul>
  )
}

export default function TermSearchPage({
  searchParams
}: {
  searchParams: {
    q: string | string[] | undefined
  },
}) {
  const router = useRouter()
  const terms = React.useMemo(() =>
    ensureArray(searchParams.q).flatMap(el => el.split(/\s+/g)).filter(el => el.length > 0),
    [searchParams.q])
  const [rawTerms, setRawTerms] = React.useState(terms.join(' '))
  const [searchTerms, setSearchTerms] = React.useState<string[]>(terms)

  return (
    <>
    <div className='flex-col'>
      <form
        className="flex flex-row items-center gap-2"
        onSubmit={evt => {
          evt.preventDefault()
          router.push(`/term-search?q=${encodeURIComponent(rawTerms)}`, {
            scroll: false,
          })
          setSearchTerms(rawTerms.split(' ').filter(el => el.length > 0))
        }}
      >
        <span className="label-text text-lg">Term</span>
        <input
          type="text"
          className="input input-bordered"
          placeholder="neuron"
          value={rawTerms}
          onChange={evt => {
            setRawTerms(evt.currentTarget.value)
          }}
          onSubmit={evt => {
            setSearchTerms(terms)
          }}
        />
        <button
          type="submit"
          className="btn normal-case"
        >Search gene sets</button>
       <div className='ml-10'> Query extracted gene set table titles to find relevant gene sets.</div>
      </form>
      <p className="prose p-2">
        try an example:
        <a
          className="font-bold text-sm cursor-pointer"
          onClick={() => {
            setRawTerms('neuron')
            setSearchTerms(['neuron'])
          }}
        > neuron</a>,
        <a
          className="font-bold text-sm cursor-pointer"
          onClick={() => {
            setRawTerms('CRISPR')
            setSearchTerms(['CRISPR'])
          }}
        > CRISPR</a>, 
        <a
          className="font-bold text-sm cursor-pointer"
          onClick={() => {
            setRawTerms('PBMC')
            setSearchTerms(['PBMC'])
          }}
        > PBMC</a>
      </p>
      </div>
      {searchTerms.length > 0 ?
      <React.Suspense fallback={<Loading/>}>
        <TermSearchResults terms={searchTerms} />
      </React.Suspense>: <></>}
    </>
  )
}
