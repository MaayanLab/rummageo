
import React from 'react'
import Image from 'next/image'
import Stats from '@/app/stats'

export default function Loading({show_total_gene_sets}: {show_total_gene_sets?: boolean}) {
  return (
    <>
        <div className="text-center p-5">
        <Image className={'rounded mx-auto'} src={'/images/loading.gif'} width={250} height={250} alt={'Loading...'}/> 
        {show_total_gene_sets ? <p>Rummaging through <Stats bold show_total_gene_sets /> that match your query.</p> : <></>}
        </div>
    </> 
  )
}




