// Deploy copy only: provenance strings in the generated data contain absolute paths of the build machine
// (/Users/<name>/Desktop/<project>/...). Rewrite them to project-relative paths in dist/families; source data is untouched.
import { readdirSync, readFileSync, statSync, writeFileSync } from 'node:fs'
import { join } from 'node:path'

const ROOT = process.argv[2] ?? 'dist/families'
const ABS = /\/(?:Users|home)\/[^/"\s]+\/(?:Desktop\/)?/g
let files = 0, hits = 0
const walk = (d) => { for (const f of readdirSync(d)) { const p = join(d, f); if (statSync(p).isDirectory()) walk(p); else if (/\.(json|csv|md|txt|sto|scores|fasta)$/.test(f)) fix(p) } }
const fix = (p) => { const s = readFileSync(p, 'utf8'); const n = (s.match(ABS) ?? []).length; if (n) { writeFileSync(p, s.replace(ABS, '')); files++; hits += n } }
walk(ROOT)
console.log(`sanitize-dist-paths: ${hits} absolute paths made project-relative in ${files} files under ${ROOT}`)
