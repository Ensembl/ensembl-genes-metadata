import {
  Accordion,
  AccordionContent,
  AccordionItem,
  AccordionTrigger,
} from "@/components/ui/accordion";

const faq = [
  {
    question: "What is Genebuild Metadata and what can I use it for?",
    answer:
      "Genebuild Metadata is a web app for exploring genome assemblies and gene annotations used by the Ensembl genebuild team and related biodiversity projects. You can quickly search assemblies and annotations by project, accession, BioProject, or taxon ID; check which assemblies are annotated or ready for annotation; and generate publication‑ready tables and figures summarising project progress, data availability, and quality.",
  },
  {
    question: "Which projects and genomes are included in the app?",
    answer:
      "The app covers all INSDC eukaryotic genome assemblies including but not limited to major biodiversity and reference genome initiatives, such as DToL, ERGA, ERGA‑pilot, EBP, VGP, HPRC, CBP, ASG, ATLASea, LACA, Rodent2K, and others. The underlying data are drawn from Ensembl genebuild tracking, ENA and related resources, and are kept in sync with the current Genebuild pipelines.",
  },
  {
    question: "What counts as a “Biodiversity project” in this app?",
    answer:
      "A biodiversity project in this app is a coordinated genome sequencing effort (for example DToL, ERGA, EBP, VGP, HPRC) that has agreed metadata definitions and project tracking within Genebuild. For these projects, the app groups assemblies and annotations by project, highlights high‑priority targets, and shows how many assemblies are unannotated, in progress, or live.",
  },
  {
    question:
      "How do I generate a report on available assemblies for a specific biodiversity project?",
    answer:
      "Go to the Report page and choose the Assemblies option. Select the biodiversity project you are interested in, and then apply any additional filters you need (for example year, assembly status, or custom groups). The app will generate summary plots and downloadable tables listing all matching assemblies, including accession IDs, taxonomy, release dates, and qualification status.",
  },
  {
    question: "How do I search for assemblies to annotate?",
    answer:
      "Open the Assemblies page. You can filter by project name, assembly accession (GCA), BioProject ID, or taxon ID. Use the “Only show current assemblies” option to hide deprecated accessions, and “Only show non‑annotated assemblies” to focus on genomes that do not yet have a Genebuild annotation. Once you click “Get Assemblies”, you will see a table with links to assembly metrics, ENA, and transcriptomic resources. Click on Assembly metrics then Get annotation candidates to further filter you assemblies.",
  },
  {
    question:
      "What does Check ENA and Check transcriptomic registry do and what will I see there?",
    answer:
      "For each assembly, the “Check ENA” toggle checks for RNAseq data submitted to ENA based on the lowest and genus taxon id. The API is reading a cache first that updates every 3 months so if your work is critical or you know there is recent new data you may need to check manually. The “Check transcriptomic registry” link takes you to the registry of transcriptomic datasets linked to that assembly or taxon, helping you confirm whether suitable RNA‑seq data are available to support annotation.",
  },
  {
    question: "How do I search for annotations?",
    answer:
      "Use the Annotations page. You can search by project name, assembly accession (GCA), BioProject ID, or taxon ID, and you can filter by annotation date. After you click “Get Annotations”, the app shows the relevant annotation records, including which assemblies have live annotations, which are in progress, and when they were generated.",
  },
  {
    question:
      "How do I request a new feature or what happens if I spot an issue?",
    answer:
      "If something looks wrong in the data or you need a new feature, please contact the Ensembl Genebuild team through the usual internal channels. Include the URL, project or assembly accession, and a short description of the issue or request so we can investigate quickly.",
  },
];


const FAQ = () => {
  return (
    <div className="px-6 py-20">
      <div className="mx-auto w-full max-w-xl">
        <h2 className="font-medium text-4xl leading-[1.15]! tracking-[-0.04em] md:text-[2.75rem]">
          Questions & Answers
        </h2>

        <Accordion
          className="mt-6"
          defaultValue={["question-0"]}
          type="multiple"
        >
          {faq.map(({ question, answer }, index) => (
            <AccordionItem key={question} value={`question-${index}`}>
              <AccordionTrigger className="text-left text-lg">
                {question}
              </AccordionTrigger>
              <AccordionContent className="text-base text-muted-foreground">
                {answer}
              </AccordionContent>
            </AccordionItem>
          ))}
        </Accordion>
      </div>
    </div>
  );
};

export default FAQ;
