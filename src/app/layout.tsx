import './globals.css'
import React from 'react'
import type { Metadata } from 'next'
import Link from "next/link"
import { ApolloWrapper } from '@/lib/apollo/provider'
import Nav from './nav'
import Stats from './stats'
import Image from 'next/image'
import { RuntimeConfig } from '@/app/runtimeConfig'
import { Open_Sans, Roboto_Mono } from 'next/font/google'

import dynamic from "next/dynamic";
import Analytics from './analytics'

const ConsentCookie = dynamic(() => import("@/components/consentCookie"), {
  ssr: false,
});

const openSans = Open_Sans({
  subsets: ['latin'],
  display: 'swap',
  //👇 Add variable to our object
  variable: '--font-opensans',
})

//👇 Configure the object for our second font
const robotoMono = Roboto_Mono({
  subsets: ['latin'],
  display: 'swap',
  variable: '--font-roboto-mono',
})

export const metadata: Metadata = {
  title: 'RummaGEO',
  description: 'Search through automatically generated signatures from GEO',
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode,
}) {
  return (
    <html lang="en" style={{ minWidth: '580px' }} className={`${openSans.variable} ${robotoMono.variable} font-sans`}>
      <ApolloWrapper>
        <RuntimeConfig>
          <body className="min-h-screen flex flex-col">
            <header>
              <div className="navbar block text-center">
                <div className="navbar-center">
                  <ul className="menu menu-horizontal gap-3 flex text-lg">
                    <Nav />
                  </ul>
                </div>
                <div className="navbar-center ml-5">
                  <React.Suspense fallback={<span className="loading loading-ring loading-lg"></span>}>
                    <Stats bold show_sets_analyzed />
                  </React.Suspense>
                </div>
              </div>
            </header>
            <main className="flex-1 flex flex-col justify-stretch mx-8 md:mx-32">
              <React.Suspense fallback={<span className="loading loading-ring loading-lg"></span>}>
                {children}
              </React.Suspense>
            </main>
            <footer className="flex-none footer p-5 mt-5 bg-neutral-900  bg-opacity-40 text-neutral-content flex place-content-evenly">
              <div className="text-center pt-6">
                <ul>
                  <li><Link href="mailto:avi.maayan@mssm.edu" target="_blank">Contact Us</Link></li>
                  <li>
                    <Link href="https://github.com/MaayanLab/rummageo" target="_blank" rel="noopener noreferrer">
                      Source Code
                    </Link>
                  </li>
                </ul>
              </div>
              <div className="text-center">
              <p>
                  <Link href="https://labs.icahn.mssm.edu/" target="_blank" rel="noopener noreferrer">
                    <Image src={'/images/ismms_white.png'} width={150} height={250} alt={'Ma&apos;ayan Lab'}/>
                  </Link>
                </p>
              </div>
              <div className="text-center pt-5">
              <p>
                <Link href="https://labs.icahn.mssm.edu/maayanlab/" target="_blank" rel="noopener noreferrer">
                  <Image className={'rounded'} src={'/images/maayanlab_white.png'} width={125} height={250} alt={'Ma&apos;ayan Lab'}/>
                </Link>
                </p>
              </div>
            </footer>
            <Analytics />
          </body>
        </RuntimeConfig>
      </ApolloWrapper>
    </html>
  )
}
