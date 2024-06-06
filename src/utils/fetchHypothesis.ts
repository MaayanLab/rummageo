'use server'

const systemPrompt = `You are an AI hypothesis generator for RummaGEO (automatically generated signatures from GEO). You should act as a biologist in hypothesizing why a high overlap may exist between the users gene set (which they will provide a description of) and the automatically generated gene set from GEO (Gene expression omnibus).`


export default async function fetchHypothesis(userDesc: string, gseSummary: string, term: string, enrichedTerms: Record<string, string[]>, enrichedStats: Record<string, string[]>) {

    var enrichedTermsString = "";
    Object.keys(enrichedTerms).forEach((library) => { 
        enrichedTermsString += `${library}: ${enrichedTerms[library].join(', ')}\n`
    })

    const prompt = `Here are two gene sets that highly overlap. The first is from a user-submmited gene set. \
    The second is a gene set automatically computed between two conditions in a study from the Gene Expression Omnibus (GEO). \
    Based upon the term name (formatted as condition 1 vs. condition 2) and the abstract of the GEO gene set, and the user submitted description of their gene set, please hypothesize about why these two gene sets have a significant high overlap.
    You should mention both the abstract of the GEO gene set and the user submitted description of their gene set in your hypothesis. You will also be provided with enrichment results from the Enrichr database to help you generate your hypothesis which shows signfigantly overlapping functional terms from the overlapping genes of the two sets.
    For each enrichment term that appears in your response, the term should appear in the exact form it was given to you (do not exclude any words or characters from a term. For example, 
    Complement And Coagulation Cascades WP558 should appear as Complement And Coagulation Cascades WP558, not Complement And Coagulation Cascades). Also, please don't use quotes around the enriched term names.
    Gene set term 1 (from GEO): ${term}
    "up" or "dn" in this term name indicates if the genes were upregulated or downregulated in the signature
    abstract of paper for gene set term 1: ${gseSummary}

    Gene set term 2: user submitted gene set
    abstract of paper for gene set term 2: ${userDesc}
    
    Enriched Terms from overlapping genes of the two sets: 
    ${enrichedTermsString}`

    const tagLine = await fetch(`https://api.openai.com/v1/chat/completions`, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'Authorization': `Bearer ${process.env.OPENAI_API_KEY}`,
        },
        body: JSON.stringify({
            model: 'gpt-4o',
            messages: [
                { "role": "system", "content": systemPrompt },
                { "role": "user", "content": prompt }
            ],
            max_tokens: 1000,
            temperature: 0
        }),
    })


    const tagLineJson = await tagLine.json()

    var hypothesis: string = tagLineJson.choices[0].message.content
    Object.keys(enrichedStats).forEach((word, index) =>{
        if (hypothesis.includes(word)) {
            hypothesis = hypothesis.replaceAll(word, `<div className="tooltip underline italic flex-wrap inline z-50" 
            data-html="true" data-tip="Library: ${enrichedStats[word][9]}&#013;&#010;
            Rank: ${enrichedStats[word][0]}&#013;&#010;
            P-value: ${Number(enrichedStats[word][2]).toExponential(2)}&#013;&#010;
            Odds Ratio: ${Number(enrichedStats[word][3]).toFixed(4)}
            ">${word}</div>`)
        }
    })

    return hypothesis
}