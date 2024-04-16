import getItem from "../../item"
import { redirect } from 'next/navigation'

export async function GET(request: Request, { params }: { params: { id: string } }) {
  const geneSet = await getItem(params.id)
  if (!geneSet.data.geneSet) return new Response(JSON.stringify({error: 'Not Found'}), { status: 404 })
  const req = await fetch('https://g2sg.cfde.cloud/api/addGeneset', {
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/json',
    },
    method: 'POST',
    body: JSON.stringify({
      term: geneSet.data.geneSet.term,
      genes: geneSet.data.geneSet.genes.nodes.map(gene => gene.symbol),
      description: `Rummagene ${geneSet.data.geneSet.description}`,
    }),
  })
  const { text: session_id } = await req.json()
  if (!session_id) return new Response(JSON.stringify({error: 'Failed to Register Gene Set'}), { status: 500 })
  redirect(`https://g2sg.cfde.cloud/analyze/${session_id}`)
}
