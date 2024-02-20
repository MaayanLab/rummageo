'use client'
import Analytics from '@/app/analytics'
import { getCookie, setCookie } from 'typescript-cookie'

export const ConsentCookie = () => {
    if (getCookie('consentCookie') === undefined) {
        return(
            <div id="cookieConsent" className='text-center fixed bg-slate-800 rounded-md p-5 bottom-0 w-full '>
                <h1 className=' font-light text-xl'>Cookie Policy</h1>
                <div>
                       <p>Is it okay to utilize Google Analytics while visiting this website to help improve user experience?</p>
                        <button  className="border p-2 m-5 rounded-lg" onClick={() => {
                            setCookie('consentCookie', 'allow')
                            window.location.reload()
                        }}>Agree</button>
                        <button  className="border p-2 m-5 rounded-lg" onClick={() => {
                            setCookie('consentCookie', 'deny')
                            window.location.reload()
                        }}>Decline</button>
                </div>
            </div>
        )
    } else if (getCookie('consentCookie') === 'deny') {
        return <></>
    } else {
        return <Analytics />
    }
}

export default ConsentCookie