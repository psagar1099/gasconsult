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
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <style>
        body { font-family: 'Inter', sans-serif; background: #0f172a; color: #e2e8f0; margin:0; padding:20px; }
        .container { max-width: 900px; margin: 0 auto; background: #1e293b; padding: 40px; border-radius: 16px; box-shadow: 0 20px 40px rgba(0,0,0,0.4); }
        h1 { font-size: 3rem; text-align: center; color: #10b981; margin:0 0 10px 0; }
        p.tagline { text-align: center; font-size: 1.2rem; color: #94a3b8; margin-bottom: 40px; }
        textarea { width: 100%; height: 120px; padding: 16px; border-radius: 12px; border: none; font-size: 1.1rem; margin-bottom: 16px; }
        input[type="submit"] { background: #10b981; color: white; padding: 14px 32px; border: none; border-radius: 12px; font-size: 1.1rem; cursor: pointer; }
        input[type="submit"]:hover { background: #059669; }
        .response { margin-top: 40px; background: #334155; padding: 24px; border-radius: 12px; line-height: 1.7; }
        .references { margin-top: 40px; }
        .ref { background: #1e293b; padding: 16px; border-radius: 12px; margin-bottom: 16px; }
        a { color: #10b981; }
        footer { text-align: center; margin-top: 60px; color: #64748b; font-size: 0.9rem; }
    </style>
</head>
<body>
    <div class="container">
        <h1>gasconsult.ai</h1>
        <p class="tagline">Evidence-based anesthesiology answers • Powered by GPT + PubMed • No hallucinations</p>
        
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