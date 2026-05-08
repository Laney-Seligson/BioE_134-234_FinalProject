### Laney:
<div style="margin-left: 20px;">
<b> 1. create_construction_file.py:</b> 
create_construction_file.py has two separate abilities:
<div style="margin-left: 20px;">
a) Can create a full construction file for a user who wants one for their specific inputs, and can take information from prior conversation in the same chat.

  Takes information already provided by the user in the conversation or prompts the user for missing information. Generates a construction file easily viewable inline to the user in the following format:

```
PCR         primerF         primerR         insert_template        insert_pcr
PCR         vectorF         vectorR         backbone_plasmid       vector_pcr
assembly_method   vector_pcr      insert_pcr      enzyme_or_reagent      assembled_plasmid
Transform   assembled_plasmid   competent_cells   antibiotic   temperature   final_construct

dsdna       insert_template    INSERT_DNA_SEQUENCE
plasmid     backbone_plasmid   BACKBONE_SEQUENCE
oligo       primerF        SEQUENCE_FORWARD_PRIMER
oligo       primerR        SEQUENCE_REVERSE_PRIMER
oligo       vectorF        SEQUENCE_VECTOR_FORWARD
oligo       vectorR        SEQUENCE_VECTOR_REVERSE
```

The user can provide any information they want, except that, due to API model limitations, backbone sequences are too long for the model to process. Therefore, we have curated a list of plasmids that the user can ask about, and the model will provide information on the ones they may select. The user then just needs to mention the plasmid they would like to use.

Example:
>Note: Break up the entire prompt into multiple prompts in order not to overwhelm the API tokens. Also, it may call errors as it realizes you did not give it all the needed information in one single prompt, but the API should continue to prompt you for the missing information.
```
You: I want to make a construction file. I want to do Golden Gate, my construct name is pET28a_REP24, my insert name is REP24, my insert sequence is atgaaaaatgttttaatggttactacttctcatgatgttatgggtaattctaatgaaaaaactggtttatggttatctgaattaactcatccttattattctattattgataaaaatattaatattgatattgtttctattatgggtggtgaaattcctattgatcctaattctgttgctcaagaagattattataatgataaatttttagctgatgataatttaaaaaatattatgaaaaattctacttctttacgtgatgttaatattaaagaatatgatgctattatttttgctggtggtcatggtactatgtgggattttcctaataatgctaatattcattctaaagttttagatatttatgctaaaaatggtgttattggtgctatttgtcatggtgttgctgctttaattaatgttaaagataataatggtcaaaatattattcgtgataaagaagttactggtttttctaataatgaagaaaaaattgttggtttaactgatgttgttcctttttctttagaagattctttagttgaagctggtgctaaatattcttctgcttctgaatggcaatcttatgttaaatctgattctaaaattattactgctcaaaatcctcaatctgctactgattttgctaaagctattaaacaatctttatttaat
Gemini: Should ask for the missing information
You: my backbone name is pET28a, and I want to use the plasmid pET28a for the backbone sequence, my insert forward primer name is repF, my insert forward primer sequence is ccataGGTCTCaATGAAAAATGTTTTAATGGTTACTA, my insert reverse primer name is repR, my insert reverse primer sequence is cagatGGTCTCaCGAGATTAAATAAAGATTGTTTAAT, my vector forward primer name is vecF
Gemini: Should ask for the missing information
You: my vector reverse primer name is vecR, my vector forward primer sequence is ccataGGTCTCaCTCGAGCACCACCACCACCACCACT, my vector reverse primer sequence is cagatGGTCTCaTCATGGTATATCTCCTTCTTAAAGT, my restriction enzyme is BsaI. Output my construction file
Gemini: Should output the construction file
```

b) Create a shorthand construction file from a paper
  The user can inquire about more information about what papers are available to choose from. They can then ask for a shorthand construction file to be created about that paper. This tool is then able to generate a shorthand construction file with the paper’s details.

Example:
```
You: What papers are available for me to look at?
Gemini: list of available papers. May ask you if you want to create a shorthand construction file
You: Yes, I want a shorthand construction file for Miao 2013
Gemini: Should output the shorthand construction file
```
</div>

<b>2. Create_construction_file.json:</b>
  This file is the C9 JSON wrapper for the construction file generator. It defines how the MCP framework and Gemini should understand and call the construction tool, including its name, description, inputs, outputs, examples, and execution details. The wrapper tells the model that this tool can generate a full sequence-based construction file, create a structured paper-information record, or generate a shorthand workflow from paper-derived metadata, depending on the selected input mode.

<b>3. validate_construction_file.py:</b>
This file contains the Python implementation of the construction file validation tool. Its job is to check whether a proposed cloning workflow is biologically consistent. It validates things like whether primers anneal correctly to the intended templates, whether they are oriented properly, whether a plausible amplicon can be formed, and whether the overall workflow structure is biologically reasonable. It can retrieve information from the construction file previously generated. This tool is useful after generating a construction file because it provides a second layer of error checking before the workflow is trusted for downstream use.

Supported step types and what is validated:
- **PCR**: primer annealing, orientation, and predicted amplicon length
- **GoldenGate**: BsaI overhang compatibility for circular assembly
- **Gibson**: overlap length and sequence compatibility at both junctions

Example prompt (happy path):
>Note: Have the API create a construction file first. For this example, let’s assume there are 2 PCR steps and a Gibson step in the construction file. In this case, the construction file is valid.
```
You: Validate this
Gemini: PCR step 1 passed, PCR step 2 passed, Gibson passed. Overall: Pass
```
Example prompt (sad path):
>Note: Have the API create a construction file first. For this example, let’s assume there are 2 PCR steps and a goldengate step in the construction file. In this case, the goldengate overhangs do not overlap.
```
You: Validate this
Gemini: PCR step 1 passed, PCR step 2 passed, goldengate failed due to not aligned overhangs. Overall: Fail
```

<b>4. validate_construction_file.json:</b>
This file is the JSON wrapper for the construction file validator. It describes the validator to the MCP framework by specifying the tool metadata, accepted inputs, expected outputs, and how the Python implementation should be executed. Its main purpose is to let Gemini recognize when the user is asking to check, verify, or debug a cloning design and then correctly route that request to the validation logic.

<b>5. Get_paper_info.py:</b>
This file contains the Python implementation of the paper information loader. Instead of forcing the model to guess details from a paper title alone, this tool reads a curated JSON record from the module’s data directory and returns structured metadata about a paper, such as organism, system, vectors, enzymes, assembly method, delivery method, validation methods, and major constraints. This makes the paper-based workflow much more reliable because downstream tools can use structured paper information rather than hallucinating experimental details.

<b>6. Get_paper_info.json:</b>
This file is the JSON wrapper for the paper information loader. It defines how Gemini can call the tool by providing a paper_id, and it tells the framework where the corresponding Python implementation lives. The purpose of this wrapper is to expose curated paper metadata as a callable MCP tool so that other tools, such as shorthand workflow generation, can use literature-derived information in a structured and reproducible way.

<b>Citations:</b>
  - Hall, Bradford, Andrew Cho, Advait Limaye, Kyoungin Cho, Jaspal Khillan, and Ashok B. Kulkarni. 2018. “Genome Editing in Mice Using CRISPR/Cas9 Technology.” Current Protocols in Cell Biology 81 (1): e57. https://doi.org/10.1002/cpcb.57.
  - Miao, Jin, Dongshu Guo, Jinzhe Zhang, Qingpei Huang, Genji Qin, Xin Zhang, Jianmin Wan, Hongya Gu, and Li-Jia Qu. 2013. “Targeted Mutagenesis in Rice Using CRISPR-Cas System.” Cell Research 23 (10): 1233–36. https://doi.org/10.1038/cr.2013.123.
  - Hall, Bradford, Andrew Cho, Advait Limaye, Kyoungin Cho, Jaspal Khillan, and Ashok B. Kulkarni. 2018. “Genome Editing in Mice Using CRISPR/Cas9 Technology.” Current Protocols in Cell Biology 81 (1): e57. https://doi.org/10.1002/cpcb.57.
  - “Super Simulator - SynBio Project Tutorials.” n.d. https://ucb-bioe-anderson-lab.github.io/cloning-tutorials/tools/supersimulator/#step-1-paste-your-construction-file.
  - New England Biolabs. n.d. “Golden Gate Assembly Domestication Tutorial.” https://www.neb.com/en-us/applications/cloning-and-synthetic-biology/dna-assembly-and-cloning/golden-gate-assembly?srsltid=AfmBOopg3N59s77LuEIArVnb6sf9Mxs3QQIvBPdP0-s_u6WLgFnhUjxh.
  - New England Biolabs. n.d. “Introduction to Gibson AssemblyTM.” https://www.neb.com/en-us/applications/cloning-and-synthetic-biology/dna-assembly-and-cloning/gibson-assembly?srsltid=AfmBOoojQV5Vac9DwTV5iHmNu7VOOlzIHbsOp5k2wHydVIZj5Fqrerij.
</div>
