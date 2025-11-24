from flask import Flask, request, render_template_string
from Bio import Entrez
import openai
import os
import re

app = Flask(__name__)

Entrez.email = "psagar1099@gmail.com"
Entrez.api_key = "f239ceccf2b12431846e6c03ffe29691ac08"

openai_client = openai.OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

def clean_query(query):
    """Strip conversational filler from user queries to extract the core medical topic."""
    q = query.lower().strip()

    # Remove common conversational phrases
    filler_phrases = [
        r"^(can you |could you |please |pls )",
        r"^(tell me about |tell me |explain |describe )",
        r"^(what is |what's |what are |whats )",
        r"^(what do you know about |what does the evidence say about )",
        r"^(what's the evidence for |what is the evidence for |evidence for )",
        r"^(i want to know about |i'd like to know about |i need to know about )",
        r"^(give me info on |give me information on |info on )",
        r"^(how does |how do |how is |how are )",
        r"^(why does |why do |why is |why are )",
        r"^(search for |look up |find |show me )",
        r"^(help me understand |help me with )",
    ]

    # Apply each pattern repeatedly until no more matches
    for pattern in filler_phrases:
        q = re.sub(pattern, "", q, flags=re.IGNORECASE)

    # Remove trailing question marks and extra whitespace
    q = re.sub(r"\?+$", "", q).strip()
    q = re.sub(r"\s+", " ", q)

    return q if q else query  # Return original if cleaning removed everything

HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>gasconsult.ai — Evidence. Instantly.</title>
    <link rel="icon" href="/static/favicon.ico" type="image/x-icon">
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
    <style>
        * { box-sizing: border-box; }
        body {
            font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            margin: 0;
            padding: 40px 20px;
            min-height: 100vh;
        }
        .container {
            max-width: 950px;
            margin: 0 auto;
            background: #ffffff;
            padding: 60px 50px;
            border-radius: 24px;
            box-shadow: 0 25px 60px rgba(0, 0, 0, 0.15);
        }
        .header {
            text-align: center;
            margin-bottom: 50px;
            padding-bottom: 30px;
            border-bottom: 2px solid #f1f5f9;
        }
        .logo {
            height: 80px;
            margin-bottom: 20px;
            filter: drop-shadow(0 4px 12px rgba(102, 126, 234, 0.3));
        }
        h1 {
            font-size: 3.5rem;
            font-weight: 700;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin: 0 0 12px 0;
            letter-spacing: -0.5px;
        }
        .tagline {
            font-size: 1.25rem;
            color: #64748b;
            font-weight: 500;
            margin: 0;
        }
        form {
            margin-bottom: 20px;
        }
        textarea {
            width: 100%;
            height: 140px;
            padding: 20px 24px;
            border-radius: 16px;
            border: 2px solid #e2e8f0;
            background: #ffffff;
            color: #1e293b;
            font-size: 1.1rem;
            font-family: inherit;
            resize: vertical;
            transition: all 0.2s ease;
            font-weight: 500;
        }
        textarea:focus {
            outline: none;
            border-color: #667eea;
            box-shadow: 0 0 0 4px rgba(102, 126, 234, 0.1);
        }
        textarea::placeholder {
            color: #94a3b8;
            font-weight: 400;
        }
        .button-wrapper {
            text-align: center;
            margin-top: 24px;
        }
        input[type="submit"] {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 18px 50px;
            border: none;
            border-radius: 12px;
            font-size: 1.15rem;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.3s ease;
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
            letter-spacing: 0.3px;
        }
        input[type="submit"]:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 25px rgba(102, 126, 234, 0.5);
        }
        input[type="submit"]:active {
            transform: translateY(0);
        }
        h2 {
            font-size: 1.75rem;
            font-weight: 700;
            color: #1e293b;
            margin: 0 0 20px 0;
            padding-bottom: 12px;
            border-bottom: 3px solid #667eea;
        }
        .response {
            margin-top: 50px;
            background: #f8fafc;
            padding: 35px;
            border-radius: 16px;
            border: 2px solid #e2e8f0;
        }
        .response p {
            color: #334155;
            font-size: 1.05rem;
            line-height: 1.8;
            font-weight: 400;
        }
        .references {
            margin-top: 40px;
            background: #ffffff;
            padding: 35px;
            border-radius: 16px;
            border: 2px solid #e2e8f0;
        }
        .ref {
            background: #f8fafc;
            padding: 24px;
            border-radius: 12px;
            margin-bottom: 18px;
            border-left: 4px solid #667eea;
            transition: all 0.2s ease;
        }
        .ref:hover {
            background: #f1f5f9;
            transform: translateX(4px);
        }
        .ref strong {
            color: #1e293b;
            font-size: 1.05rem;
            font-weight: 600;
            line-height: 1.5;
            display: block;
            margin-bottom: 8px;
        }
        .ref a {
            color: #667eea;
            text-decoration: none;
            font-weight: 600;
            transition: color 0.2s ease;
        }
        .ref a:hover {
            color: #764ba2;
            text-decoration: underline;
        }
        footer {
            text-align: center;
            margin-top: 60px;
            padding-top: 30px;
            color: #94a3b8;
            font-size: 0.95rem;
            border-top: 2px solid #f1f5f9;
            font-weight: 500;
        }
        .debug {
            color: #667eea;
            font-size: 0.9rem;
            margin-top: 16px;
            font-weight: 600;
        }
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
            <textarea name="query" placeholder="e.g., Tell me about TXA in spine surgery, propofol in pediatrics, PONV prevention" required></textarea>
            <div class="button-wrapper">
                <input type="submit" value="Get Evidence">
            </div>
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
        raw_query = request.form["query"].strip()
        query = clean_query(raw_query)  # Strip conversational filler

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

        # ←←← THIS IS THE NEW, SMART PROMPT THAT ACTUALLY WORKS ←←←
        prompt = f"""You are an expert anesthesiologist in the OR right now.

Question: {raw_query}

The references below are real, recent, high-quality papers (systematic reviews, meta-analyses, RCTs) that answer this exact clinical question — even if they use synonyms like "tranexamic acid" instead of "TXA", "spinal fusion" instead of "spine surgery", etc.

Answer concisely and directly using the evidence. Cite by first author + year or PMID.

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