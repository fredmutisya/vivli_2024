# vivli_2024- InsightCAre.ai

This code was developed by  Dr. Rachael Kanguha and Dr. Fredrick Mutisya using Pfizer-Atlas, Venatorx-GEARS and Paratek-Keystone data as part of the 2024 Vivli data challenge. For more information on Vivli and the challenge, visit the link 
[Vivli AMR](https://amr.vivli.org/)

![AntiMicro.ai logo](https://github.com/fredmutisya/vivli_2024/website/1.jpg)


### InsightCare.AI: Streamlit Application Guide

Welcome to the InsightCare.AI Streamlit application. This guide will walk you through the setup, functionality, and features of the InsightCare.AI platform, which aims to enhance healthcare decision-making through AI-powered data analytics and visualization.

#### Setting Up the Environment

First, ensure you have all necessary libraries installed. You can install the required libraries using the following pip command:

```bash
pip install streamlit pandas matplotlib Bio openai requests streamlit_lottie Pillow
```

Next, set up your environment variables. You will need an API key for OpenAI and your email for the NCBI Entrez API:

1. Create a `.env` file in your project directory.
2. Add the following lines to your `.env` file, replacing `your_openai_api_key` and `your_email@example.com` with your actual OpenAI API key and email:

```plaintext
OPENAI_API_KEY=your_openai_api_key
ENTREZ_EMAIL=your_email@example.com
```

#### Running the Application

To start the application, navigate to your project directory in the terminal and run:

```bash
streamlit run Hello.py
```


#### Application Structure

The InsightCare.AI application is structured into several key areas:

1. **Clinical Surveillance**
   - View and interact with clinical data such as patient histories, diagnoses, and treatments.
   - Generate dynamic reports by filtering data based on demographics, conditions, or treatment outcomes.
   - Utilize AI to predict trends and possible outbreaks based on current data entries.

2. **Laboratory Surveillance**
   - Analyze antibiogram data to understand resistance patterns across various pathogens and antibiotics.
   - Visualize data through interactive charts and graphs to better interpret the spread and resistance levels of different bacteria.
   - Use predictive models to forecast resistance trends and advise on antibiotic stewardship programs.

3. **Consult Research**
   - Access PubMed directly from the application to search for relevant medical studies using dynamically generated queries based on user input.
   - Summarize research articles using AI, providing quick insights and relevant information to users without the need for deep reading.
   - Generate citations and bibliographies automatically for research materials to aid in report preparation and academic writing.

Each of these areas is accessible through tabs within the Streamlit application.

#### Key Features

- **Dynamic Background**: The application features a grey-to-blue gradient background that changes subtly.
- **Interactive Sidebar**: Includes patient lists and Lottie animations for better user engagement.
- **Data Visualization**: Use Matplotlib for generating plots and Streamlit components for displaying interactive elements.

#### Sample Code Snippets

Here are a few snippets from the application's codebase:

**Load CSV and Set Up Entrez Email**:
```python
import pandas as pd
from Bio import Entrez

# Load data
df_antibiotics = pd.read_csv('cleaned_antibiotics.csv')
Entrez.email = "your_email@example.com"  # Replace with your email
```

**Search PubMed and Display Results**:
```python
from Bio import Entrez, Medline

def search_pubmed(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_details(id_list):
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
    records = list(Medline.parse(handle))
    handle.close()
    return records

# Example usage
id_list = search_pubmed("diabetes AND metformin")
details = fetch_details(id_list)
```

#### Deployment

To deploy the application, you can use Streamlit sharing, Heroku, or any other platform that supports Python apps. Ensure that you set up environment variables securely on the deployment platform.

#### Conclusion

The InsightCare.AI application leverages Streamlit, OpenAI, and NCBI's PubMed to provide a comprehensive tool for healthcare professionals. By integrating AI and advanced data analytics, it enhances clinical decision support and facilitates access to the latest research.
