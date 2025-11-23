from flask import Flask, request, render_template_string
from Bio import Entrez
import openai
import os

app = Flask(__name__)

Entrez.email = "psagar1099@gmail.com"
Entrez.api_key = "f239ceccf2b12431846e6c03ffe29691ac08"

openai_client = openai.OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>gasconsult.ai — Evidence. Instantly.</title>
    <link rel="icon" href="/static/favicon.ico" type="image/x-icon">
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <style>
        body { font-family: 'Inter', sans-serif; background:#0f172a; color:#e2e8f0; margin:0; padding:20px 0; }
        .container { max-width:900px; margin:0 auto; background:#1e293b; padding:40px; border-radius:16px; box-shadow:0 20px 40px rgba(0,0,0,0.4); }
        .header { text-align:center; margin-bottom:40px; }
        .logo { height:90px; filter:drop-shadow(0 0 20px rgba(16,185,129,0.4)); }
        h1 { font-size:3rem; color:#10b981; margin:10px 0; }
        .tagline { font-size:1.3rem; color:#94a3b8; margin:16px 0 40px; }
        textarea { width:100%; height:130px; padding:18px; border-radius:16px; border:none; background:#334155; color:white; font-size:1.1rem; }
        input[type="submit"] { background:#10b981; color:white; padding:16px 40px; border:none; border-radius:16px; font-size:1.2rem; cursor:pointer; }
        input[type="submit"]:hover { background:#059669; transform:translateY(-2px); }
        .response, .references { margin-top:40px; background:#1e293b; padding:28px; border-radius:16px; }
        .ref { background:#0f172a; padding:20px; border-radius:12px; margin-bottom:16px; }
        footer { text-align:center; margin-top:80px; color:#64748b; font-size:0.9rem; }
        .debug { color: #10b981; font-size: 0.9rem; margin-top: 10px; }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <img src="/static/logo.png" alt="gasconsult.ai" class="logo">
            <h1>gasconsult.ai</h1>
            <p class="tagline">Evidence-based anesthesiology answers • No hallucinations</p>
        </div>
        
        <form method="post">
            <textarea name="query" placeholder="e.g., TXA in spine surgery, blood loss scoliosis, propofol peds" required></textarea>
            <center><input type="submit" value="Get Evidence"></center>
        </form>

        {% if answer %}
        <div class="response">
            <h2>Answer</h2>
            {{ answer|safe }}
            <div class="debug">Debug: Fetched {{ num_papers }} papers</div>
        </div>
        <div class="references">
            <h2>References</h2>
            {% for ref in refs %}
            <div class="ref">
                <strong>{{ ref.title }}</strong><br>
                {{ ref.authors }} • {{ ref.journal }} ({{ ref.year }})<br>
                <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank">PubMed →</a>
            </div>
            {% endfor %}
        </div>
        {% endif %}
        
        <footer>Not medical advice • For educational use only • © 2025 gasconsult.ai</footer>
    </div>
</body>
</html>
"""

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        query = request.form["query"].strip()

        q = query.lower()
        q = q.replace("txa", '"tranexamic acid" OR TXA')
        q = q.replace("blood loss", '"blood loss" OR hemorrhage OR transfusion')
        q = q.replace("spine surgery", '"spine surgery" OR "spinal fusion" OR scoliosis')
        q = q.replace("pediatric", 'pediatric OR children OR peds')
        q = q.replace("ponv", 'PONV OR "postoperative nausea"')

        search_term = (
            f'({q}) AND '
            f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
            f'"Cochrane Database Syst Rev"[ta] OR guideline[pt]) AND '
            f'("2015/01/01"[PDAT] : "3000"[PDAT])'
        )

        # Try anesthesiology first
        handle = Entrez.esearch(db="pubmed", term=f'anesthesiology[MeSH Terms] AND {search_term}', retmax=15, sort="relevance", api_key=Entrez.api_key)
        result = Entrez.read(handle)
        ids = result["IdList"]

        if not ids:
            handle = Entrez.esearch(db="pubmed", term=search_term, retmax=15, sort="relevance", api_key=Entrez.api_key)
            result = Entrez.read(handle)
            ids = result["IdList"]

        if not ids:
            return render_template_string(HTML, answer="<p>No high-quality recent evidence found. Try rephrasing.</p>", num_papers=0, refs=[])

        handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml", api_key=Entrez.api_key)
        papers = Entrez.read(handle)["PubmedArticle"]

        refs = []
        context = ""
        for p in papers[:12]:
            try:
                art = p["MedlineCitation"]["Article"]
                title = art.get("ArticleTitle", "No title")
                abstract = " ".join(str(t) for t in art.get("Abstract", {}).get("AbstractText", [])) if art.get("Abstract") else ""
                authors = ", ".join([a.get("LastName","") + " " + (a.get("ForeName","")[:1]+"." if a.get("ForeName") else "") for a in art.get("AuthorList",[])[:5]])
                journal = art["Journal"].get("Title", "Unknown")
                year = art["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A")
                pmid = p["MedlineCitation"]["PMID"]

                refs.append({"title": title, "authors": authors, "journal": journal, "year": year, "pmid": pmid})
                context += f"Title: {title}\nAbstract: {abstract}\nAuthors: {authors}\nJournal: {journal} ({year})\nPMID: {pmid}\n\n"
            except:
                continue

        num_papers = len(refs)

        prompt = f"""You are an expert anesthesiologist. Answer the question using ONLY the references below. Cite by title or PMID. Be concise and direct.

Question: {query}

References:
{context}

Answer:"""

        try:
            response = openai_client.chat.completions.create(
                model="gpt-4o",
                messages=[{"role": "user", "content": prompt}],
                temperature=0.1
            ).choices[0].message.content
        except Exception as e:
            response = f"AI error: {str(e)}"

        return render_template_string(HTML, answer=response, refs=refs, num_papers=num_papers)

    return render_template_string(HTML)

if __name__ == "__main__":
    app.run(debug=True)