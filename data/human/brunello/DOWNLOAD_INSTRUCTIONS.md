# Instructions for LLM: Download Brunello sgRNA Library Sequences Zip

## Objective
Download the official "Brunello_library_v1.1_sequences.txt.zip" file from Addgene, which contains the sgRNA sequences for the Brunello genome-wide knockout library (Doench et al., 2016, Nature Biotechnology, DOI: 10.1038/nbt.3437).

## Step-by-Step Instructions

1. **Search for the Brunello Library Page on Addgene**:
   - Open a web browser and navigate to Addgene's search or pooled libraries section.
   - Search for "Brunello" or go directly to the known URL: https://www.addgene.org/pooled-library/broadgpp-brunello/
   - Confirm the page is for "Brunello CRISPR Knockout Pooled Library" by the Broad Institute.

2. **Locate the Sequences Download Link**:
   - On the library page, look for a section titled "Sequences" or "Download Sequences".
   - The file to download is named "Brunello_library_v1.1_sequences.txt.zip" (approximately 2.6 MB).
   - If the direct link is not visible, check the "Files" or "Supplementary Materials" tab.

3. **Handle Authentication and Terms**:
   - Addgene may require you to log in or agree to terms of use. If prompted:
     - Create a free Addgene account if you don't have one (use a valid email).
     - Accept any terms related to academic use, non-commercial distribution, and citation requirements.
   - If the download is restricted, note that and suggest contacting Addgene support for access.

4. **Download the File**:
   - Click the download link for "Brunello_library_v1.1_sequences.txt.zip".
   - Save the file to a secure location (e.g., your downloads folder).
   - Do not rename the file; keep it as "Brunello_library_v1.1_sequences.txt.zip".

5. **Verify the Download**:
   - Check the file size: It should be around 2.6 MB.
   - Attempt to unzip it (e.g., using `unzip` command or a zip tool). It should contain "Brunello_library_v1.1_sequences.txt" (a tab-separated values file with sgRNA sequences).
   - If the file is corrupted or not a zip, retry the download or check for alternative sources.

6. **Alternative Sources if Addgene Fails**:
   - If Addgene download is blocked, try the paper's supplementary materials:
     - Go to the Nature Biotechnology article: https://www.nature.com/articles/nbt.3437
     - Look for "Supplementary Information" or "Data Availability".
     - Download the zip from the journal's site (may require Springer Nature login).
   - Another option: Search academic repositories like Zenodo or Figshare for "Brunello library sequences" and download from there if available.
   - **Broader Search Strategies**:
     - Use search engines (Google, DuckDuckGo) with queries like: "Brunello_library_v1.1_sequences.txt.zip download", "Doench 2016 sgRNA sequences zip", "Brunello CRISPR library sequences".
     - Check GitHub repositories: Search GitHub for "Brunello" and look for forks or uploads of the sequences (e.g., in CRISPR-related repos).
     - Academic databases: Search PubMed, NCBI GEO, or Dryad for datasets associated with the paper DOI.
     - Community forums: Check ResearchGate, BioStars, or Reddit (e.g., r/bioinformatics) for shared links or mirrors.
     - If found on a non-official site, verify the file integrity by checking file size (~2.6 MB) and SHA256 against known values (if available online).
     - Contact authors: Email the corresponding author of the paper (Feng Zhang or John Doench) for direct access if all else fails.

7. **Post-Download Actions**:
   - Compute the SHA256 checksum of the downloaded zip file for verification (use tools like `shasum -a 256 filename.zip`).
   - Extract the contents and inspect the TSV file to ensure it has columns like Gene, sgRNA ID, and Sequence.
   - If successful, the file is ready for use in data curation scripts.

## Important Notes
- **Citation Requirement**: Always cite the original paper: Doench et al., 2016, Nature Biotechnology.
- **License**: The data is provided under Addgene's terms; check for any restrictions on redistribution.
- **Troubleshooting**: If downloads fail due to geo-restrictions or authentication, use a VPN or contact the source directly.
- **Security**: Ensure you're downloading from official sources to avoid malware.

## Expected Outcome
- A valid zip file containing the Brunello sequences TSV, ready for extraction and processing.
