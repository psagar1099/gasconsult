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
        r"(indications to use |indications for |indication for |when to use )",
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
    <title>gasconsult.ai — Evidence-Based Anesthesiology</title>
    <link rel="icon" href="/static/favicon.svg" type="image/svg+xml">
    <link rel="alternate icon" href="/static/favicon.ico">
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'SF Pro Text', 'Helvetica Neue', Arial, sans-serif;
            background: #fafafa;
            color: #1d1d1f;
            line-height: 1.6;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
        }

        /* Navigation */
        nav {
            background: rgba(255, 255, 255, 0.8);
            backdrop-filter: saturate(180%) blur(20px);
            border-bottom: 1px solid rgba(0, 0, 0, 0.1);
            padding: 20px 0;
            position: sticky;
            top: 0;
            z-index: 1000;
        }

        nav .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 0 40px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }

        .logo-text {
            font-size: 1.5rem;
            font-weight: 600;
            background: linear-gradient(135deg, #0071e3 0%, #00c3ff 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            letter-spacing: -0.5px;
        }

        .logo-symbol {
            display: inline-block;
            margin-right: 8px;
            font-size: 1.8rem;
        }

        /* Hero Section */
        .hero {
            max-width: 980px;
            margin: 80px auto 60px;
            padding: 0 40px;
            text-align: center;
        }

        h1 {
            font-size: 4rem;
            font-weight: 700;
            letter-spacing: -1.5px;
            margin-bottom: 20px;
            color: #1d1d1f;
            line-height: 1.1;
        }

        .subtitle {
            font-size: 1.5rem;
            color: #6e6e73;
            font-weight: 400;
            margin-bottom: 50px;
            line-height: 1.5;
        }

        /* Search Section */
        .search-container {
            max-width: 880px;
            margin: 0 auto 80px;
            padding: 0 40px;
        }

        .search-box {
            background: white;
            border-radius: 20px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
            padding: 40px;
            transition: all 0.3s ease;
        }

        .search-box:hover {
            box-shadow: 0 8px 30px rgba(0, 0, 0, 0.12);
        }

        textarea {
            width: 100%;
            min-height: 120px;
            padding: 20px;
            font-size: 1.1rem;
            font-family: inherit;
            border: 2px solid #e5e5e7;
            border-radius: 12px;
            resize: vertical;
            transition: all 0.2s ease;
            background: #fafafa;
            color: #1d1d1f;
        }

        textarea:focus {
            outline: none;
            border-color: #0071e3;
            background: white;
            box-shadow: 0 0 0 4px rgba(0, 113, 227, 0.1);
        }

        textarea::placeholder {
            color: #86868b;
        }

        .button-wrapper {
            margin-top: 24px;
            text-align: center;
        }

        input[type="submit"] {
            background: #0071e3;
            color: white;
            font-size: 1.1rem;
            font-weight: 500;
            padding: 16px 48px;
            border: none;
            border-radius: 980px;
            cursor: pointer;
            transition: all 0.3s ease;
            letter-spacing: -0.2px;
        }

        input[type="submit"]:hover {
            background: #0077ed;
            transform: scale(1.02);
            box-shadow: 0 4px 20px rgba(0, 113, 227, 0.3);
        }

        input[type="submit"]:active {
            transform: scale(0.98);
        }

        /* Results Section */
        .results-container {
            max-width: 980px;
            margin: 0 auto 80px;
            padding: 0 40px;
        }

        .response {
            background: white;
            border-radius: 20px;
            padding: 50px;
            margin-bottom: 40px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
        }

        h2 {
            font-size: 2rem;
            font-weight: 700;
            color: #1d1d1f;
            margin-bottom: 24px;
            letter-spacing: -0.5px;
        }

        .response p {
            font-size: 1.1rem;
            line-height: 1.8;
            color: #1d1d1f;
            margin-bottom: 16px;
        }

        .references {
            background: white;
            border-radius: 20px;
            padding: 50px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
        }

        .ref {
            padding: 28px;
            background: #f5f5f7;
            border-radius: 12px;
            margin-bottom: 20px;
            border-left: 4px solid #0071e3;
            transition: all 0.2s ease;
        }

        .ref:hover {
            background: #e8e8ed;
            transform: translateX(4px);
        }

        .ref strong {
            display: block;
            font-size: 1.1rem;
            font-weight: 600;
            color: #1d1d1f;
            margin-bottom: 12px;
            line-height: 1.5;
        }

        .ref-meta {
            color: #6e6e73;
            font-size: 0.95rem;
            margin-bottom: 10px;
            line-height: 1.6;
        }

        .ref a {
            color: #0071e3;
            text-decoration: none;
            font-weight: 500;
            transition: color 0.2s ease;
        }

        .ref a:hover {
            color: #0077ed;
            text-decoration: underline;
        }

        /* Footer */
        footer {
            max-width: 980px;
            margin: 80px auto 40px;
            padding: 40px 40px 0;
            text-align: center;
            border-top: 1px solid #d2d2d7;
            color: #86868b;
            font-size: 0.9rem;
        }

        .debug {
            display: inline-block;
            background: #f5f5f7;
            color: #6e6e73;
            padding: 8px 16px;
            border-radius: 12px;
            font-size: 0.85rem;
            margin-top: 20px;
            font-weight: 500;
        }

        /* Responsive Design */
        @media (max-width: 768px) {
            h1 {
                font-size: 2.5rem;
            }

            .subtitle {
                font-size: 1.2rem;
            }

            .hero {
                margin: 40px auto 40px;
            }

            .search-box,
            .response,
            .references {
                padding: 30px;
            }

            nav .container,
            .hero,
            .search-container,
            .results-container,
            footer {
                padding-left: 20px;
                padding-right: 20px;
            }
        }
    </style>
</head>
<body>
    <nav>
        <div class="container">
            <div class="logo-text">
                <span class="logo-symbol">⚕</span>gasconsult.ai
            </div>
        </div>
    </nav>

    <div class="hero">
        <h1>Evidence-Based<br>Anesthesiology Answers</h1>
        <p class="subtitle">Get instant access to peer-reviewed research with AI-powered synthesis. No hallucinations, just evidence.</p>
    </div>

    <div class="search-container">
        <div class="search-box">
            <form method="post">
                <textarea name="query" placeholder="Ask anything about anesthesiology... e.g., Tell me about TXA in spine surgery" required></textarea>
                <div class="button-wrapper">
                    <input type="submit" value="Get Evidence">
                </div>
            </form>
        </div>
    </div>

    {% if answer %}
    <div class="results-container">
        <div class="response">
            <h2>Answer</h2>
            {{ answer|safe }}
            <div class="debug">📊 Fetched {{ num_papers }} papers from PubMed</div>
        </div>

        <div class="references">
            <h2>References</h2>
            {% for ref in refs %}
            <div class="ref">
                <strong>{{ ref.title }}</strong>
                <div class="ref-meta">
                    {{ ref.authors }}<br>
                    {{ ref.journal }} • {{ ref.year }}
                </div>
                <a href="https://pubmed.ncbi.nlm.nih.gov/{{ ref.pmid }}/" target="_blank">View on PubMed →</a>
            </div>
            {% endfor %}
        </div>
    </div>
    {% endif %}

    <footer>
        <p>Not medical advice • For educational use only<br>
        © 2025 gasconsult.ai</p>
    </footer>
</body>
</html>
"""

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        raw_query = request.form["query"].strip()
        query = clean_query(raw_query)  # Strip conversational filler

        # Expand synonyms and medical abbreviations
        q = query.lower()
        q = q.replace(" versus ", " OR ")
        q = q.replace(" vs ", " OR ")
        q = q.replace(" vs. ", " OR ")
        q = q.replace(" over ", " OR ")
        q = q.replace(" compared to ", " OR ")
        q = q.replace(" compared with ", " OR ")
        q = q.replace("txa", '"tranexamic acid" OR TXA')
        q = q.replace("blood loss", '"blood loss" OR hemorrhage OR transfusion')
        q = q.replace("spine surgery", '"spine surgery" OR "spinal fusion" OR scoliosis')
        q = q.replace("peds", 'pediatric OR children OR peds')
        q = q.replace("pediatric", 'pediatric OR children OR peds')
        q = q.replace("ponv", 'PONV OR "postoperative nausea"')
        q = q.replace("propofol", '"propofol"[MeSH Terms] OR propofol')
        q = q.replace("etomidate", '"etomidate"[MeSH Terms] OR etomidate')

        search_term = (
            f'({q}) AND '
            f'(systematic review[pt] OR meta-analysis[pt] OR "randomized controlled trial"[pt] OR '
            f'"Cochrane Database Syst Rev"[ta] OR guideline[pt]) AND '
            f'("2015/01/01"[PDAT] : "3000"[PDAT])'
        )

        # Try anesthesiology-specific high-quality evidence first
        handle = Entrez.esearch(db="pubmed", term=f'anesthesiology[MeSH Terms] AND {search_term}', retmax=15, sort="relevance", api_key=Entrez.api_key)
        result = Entrez.read(handle)
        ids = result["IdList"]

        # Fallback 1: Try without anesthesiology restriction
        if not ids:
            handle = Entrez.esearch(db="pubmed", term=search_term, retmax=15, sort="relevance", api_key=Entrez.api_key)
            result = Entrez.read(handle)
            ids = result["IdList"]

        # Fallback 2: Drop publication type restrictions, keep recent papers
        if not ids:
            broader_search = f'({q}) AND ("2015/01/01"[PDAT] : "3000"[PDAT])'
            handle = Entrez.esearch(db="pubmed", term=broader_search, retmax=15, sort="relevance", api_key=Entrez.api_key)
            result = Entrez.read(handle)
            ids = result["IdList"]

        # Fallback 3: Even broader - any recent paper in anesthesiology
        if not ids:
            broadest_search = f'({q}) AND anesthesiology AND ("2010/01/01"[PDAT] : "3000"[PDAT])'
            handle = Entrez.esearch(db="pubmed", term=broadest_search, retmax=15, sort="relevance", api_key=Entrez.api_key)
            result = Entrez.read(handle)
            ids = result["IdList"]

        if not ids:
            return render_template_string(HTML, answer="<p>No relevant evidence found. Try rephrasing your question or using different medical terms.</p>", num_papers=0, refs=[])

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