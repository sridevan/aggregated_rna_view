// Static homepage thumbnails of each family's default (display) representative RNA chain.
// Rendering helpers: scripts/thumbnail3dmol.ts (the former Structures-tab 3Dmol code): the chain's coordinates from the PDBe ModelServer
// (label_asym_id, model 1; nothing stored), all polymer atoms of the chain (no water / ligands), 3Dmol.js cartoon with nucleotides drawn
// as base ladders, principal-axes orientation, zoomTo fit.
// Rendered at 2x in headless Chromium and downscaled to WebP; writes public/thumbnails/<ACC>.webp and manifest.json (provenance).
// Usage: npm run dev   (in another terminal)   then   npm run thumbnails [-- RF00167 RF00177]
import { chromium } from '@playwright/test'
import { writeFileSync, readFileSync } from 'node:fs'
import { fileURLToPath } from 'node:url'
import { dirname, join } from 'node:path'

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)))
const BASE = process.env.DEV_URL ?? 'http://127.0.0.1:5173'
const ACCS = process.argv.slice(2).length ? process.argv.slice(2) : ['RF00167', 'RF00162', 'RF00059', 'RF00005', 'RF00177', 'RF00011']
const OUT_PX = 240, RENDER_PX = 480, COLOR = '#5B78A8'   // slate blue (one colour for all four thumbnails)

const browser = await chromium.launch()
const page = await browser.newPage({ viewport: { width: 800, height: 800 }, deviceScaleFactor: 1 })
await page.goto(`${BASE}/`)
const manifestPath = join(ROOT, 'public/thumbnails/manifest.json')
let manifest = {}
try { manifest = JSON.parse(readFileSync(manifestPath, 'utf8')) } catch { /* new */ }
for (const acc of ACCS) {
  const r = await page.evaluate(async ({ acc, OUT_PX, RENDER_PX, COLOR }) => {
    const av = await import('/scripts/thumbnail3dmol.ts')
    const fam = await (await fetch(`/families/${acc}/family.json`)).json()
    const rep = fam.display_representative ?? fam.representative
    const url = `https://www.ebi.ac.uk/pdbe/model-server/v1/${rep.pdb_id}/atoms?label_asym_id=${encodeURIComponent(rep.label_asym_id)}&model_nums=1&encoding=cif`
    const atoms = av.parseAtomSite(await (await fetch(url)).text())
    const authChains = [...new Set(atoms.map((a) => a.chain))]
    const drawn = atoms.filter((a) => !a.het)   // polymer residues only (label_seq_id set), bases included
    const c4 = drawn.filter((a) => a.name === "C4'")
    const o = av.principal(c4.map((a) => [a.x, a.y, a.z]))
    const pdb = av.toPdb(av.drawableNucleotides(drawn), o ?   // modified nucleotides drawn with their parent base name (continuous cartoon)
      (p) => { const d = [p[0] - o.c[0], p[1] - o.c[1], p[2] - o.c[2]]; return [0, 1, 2].map((i) => o.R[i][0] * d[0] + o.R[i][1] * d[1] + o.R[i][2] * d[2]) } : undefined)
    const $3Dmol = await av.load3Dmol()
    document.body.innerHTML = ''
    const host = document.createElement('div'); host.style.cssText = `position:absolute;left:0;top:0;width:${RENDER_PX}px;height:${RENDER_PX}px`; document.body.appendChild(host)
    const v = $3Dmol.createViewer(host, { backgroundColor: 'white', antialias: true })
    v.addModel(pdb, 'pdb')
    v.setStyle({}, av.cartoon(COLOR, 'normal', c4.length > 800 ? 1.4 : 1.2))   // backbone tube + nucleotide ladders (3Dmol draws the bases when base atoms are present)
    v.zoomTo(); v.zoom(0.9); v.render()
    await new Promise((res) => setTimeout(res, 400)); v.render()
    const src = host.querySelector('canvas')
    const out = document.createElement('canvas'); out.width = out.height = OUT_PX
    const ctx = out.getContext('2d'); ctx.imageSmoothingQuality = 'high'; ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, OUT_PX, OUT_PX)
    ctx.drawImage(src, 0, 0, src.width, src.height, 0, 0, OUT_PX, OUT_PX)
    // margin check: no non-white pixel on the outer 2 px border (nothing cropped)
    const px = ctx.getImageData(0, 0, OUT_PX, OUT_PX).data
    let edge = 0
    for (let y = 0; y < OUT_PX; y++) for (let x = 0; x < OUT_PX; x++) if (x < 2 || y < 2 || x >= OUT_PX - 2 || y >= OUT_PX - 2) { const i = (y * OUT_PX + x) * 4; if (px[i] < 245 || px[i + 1] < 245 || px[i + 2] < 245) edge++ }
    return { dataUrl: out.toDataURL('image/webp', 0.9), rep: { seq_name: rep.seq_name, pdb_id: rep.pdb_id, entity_id: rep.entity_id, label_asym_id: rep.label_asym_id, auth_asym_id: rep.auth_asym_id },
      authChains, nResidues: c4.length, url, edge }
  }, { acc, OUT_PX, RENDER_PX, COLOR })
  if (r.authChains.length !== 1 || r.authChains[0] !== r.rep.auth_asym_id) throw new Error(`${acc}: ModelServer returned chains ${r.authChains} for ${r.rep.label_asym_id}, expected ${r.rep.auth_asym_id}`)
  if (r.edge) throw new Error(`${acc}: structure touches the image border (${r.edge} px) - cropped`)
  const buf = Buffer.from(r.dataUrl.split(',')[1], 'base64')
  writeFileSync(join(ROOT, `public/thumbnails/${acc}.webp`), buf)
  manifest[acc] = { file: `${acc}.webp`, ...r.rep, n_residues_drawn: r.nResidues, coordinate_source: r.url, representation: '3Dmol.js cartoon: backbone tube with nucleotide ladders (all polymer atoms), principal-axes orientation',
    colour: COLOR, size_px: OUT_PX, bytes: buf.length, generated: new Date().toISOString(), generator: 'frontend/scripts/make-thumbnails.mjs' }
  console.log(`${acc}: ${r.rep.pdb_id.toUpperCase()} ${r.rep.auth_asym_id} (label ${r.rep.label_asym_id}), ${r.nResidues} nt, ${buf.length} bytes`)
}
writeFileSync(manifestPath, JSON.stringify(manifest, null, 1) + '\n')
await browser.close()
