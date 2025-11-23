from flask import Flask, request, render_template_string
from Bio import Entrez
import openai
import os

app = Flask(__name__)

# ←←← CHANGE THIS TO YOUR REAL EMAIL ←←←
Entrez.email = "your-email@example.com"        # ← CHANGE THIS LINE

openai.api_key = os.getenv("OPENAI_API_KEY")

HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>gasconsult.ai — Evidence. Instantly.</title>
    <link rel="icon" href="/static/favicon.ico" type="image/x-icon">
    <link rel="shortcut icon" href="/static/favicon.ico" type="image/x-icon">
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <style>
        body { font-family: 'Inter', sans-serif; background: #0f172a; color: #e2e8f0; margin:0; }
        .container { max-width: 900px; margin: 0 auto; padding: 40px 20px; }
        .header { text-align: center; margin-bottom: 40px; }
        .logo { height: 90px; margin-bottom: 12px; filter: drop-shadow(0 0 20px rgba(16,185,129,0.4)); }
        h1 { font-size: 3rem; color: #10b981; margin:0; }
        .tagline { font-size: 1.3rem; color: #94a3b8; margin: 16px 0 40px; }
        textarea { width: 100%; height: 130px; padding: 18px; border-radius: 16px; border: none; font-size: 1.1rem; background:#1e293b; color:white; }
        input[type="submit"] { background: #10b981; color: white; padding: 16px 40px; border: none; border-radius: 16px; font-size: 1.2rem; cursor: pointer; transition:0.2s; }
        input[type="submit"]:hover { background: #059669; transform: translateY(-2px); }
        .response, .references { margin-top: 40px; background: #1e293b; padding: 28px; border-radius: 16px; }
        .ref { background: #0f172a; padding: 20px; border-radius: 12px; margin-bottom: 16px; }
        footer { text-align: center; margin-top: 80px; color: #64748b; font-size: 0.9rem; }
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
            <textarea name="query" placeholder="e.g., Best evidence for TXA in spine surgery?" required></textarea>
            <center><input type="submit" value="Get Evidence"></center>
        </form>

        {% if answer %}
        <div class="response">
            <h2>Answer</h2>
            {{ answer|safe }}
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
        query = request.form["query"]

        # Search only high-quality anesthesiology evidence
        search = f'anesthesiology[MeSH] AND ({query}) AND (systematic[sb] OR randomized controlled trial[pt] OR meta-analysis[pt])'
        handle = Entrez.esearch(db="pubmed", term=search, retmax=10, sort="relevance")
        result = Entrez.read(handle)
        ids = result["IdList"]

        if not ids:
            return render_template_string(HTML, answer="<p>No high-quality evidence found. Try broadening your question.</p>")

        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        papers = Entrez.read(handle)["PubmedArticle"]

        refs = []
        context = ""
        for p in papers[:10]:
            art = p["MedlineCitation"]["Article"]
            title = art.get("ArticleTitle", "No title")
            abstract = " ".join([t for t in art.get("Abstract", {}).get("AbstractText", [])]) if art.get("Abstract") else ""
            authors = ", ".join([a.get("LastName","") + " " + a.get("ForeName","")[:1] + "." for a in art.get("AuthorList",[])[:5]])
            journal = art["Journal"]["Title"]
            year = art["Journal"]["JournalIssue"]["PubDate"].get("Year", "Year?")
            pmid = p["MedlineCitation"]["PMID"]

            refs.append({"title": title, "authors": authors, "journal": journal, "year": year, "pmid": pmid})
            context += f"Title: {title}\nAbstract: {abstract}\n\n"

        prompt = f"You are an expert anesthesiologist. Answer ONLY using the references below. Cite titles in parentheses. If evidence is weak or absent, say so.\n\nQuestion: {query}\n\nReferences:\n{context}\n\nAnswer:"

        response = openai.ChatCompletion.create(
            model="gpt-4o",        # works today; change to gpt-5 when you have access
            messages=[{"role": "user", "content": prompt}],
            temperature=0.2
        ).choices[0].message.content

        return render_template_string(HTML, answer=response, refs=refs)

    return render_template_string(HTML)

if __name__ == "__main__":
    app.run(debug=True)
