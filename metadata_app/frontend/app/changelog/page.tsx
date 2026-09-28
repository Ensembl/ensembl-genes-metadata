import { Dot } from "lucide-react";
import { Separator } from "@/components/ui/separator";

const changelogs = [
  {
    version: "1.0.0",
    date: "2026-01-01",
    title: "Initial Launch",
    description:
      "We're live! The very first release of our app is here. This version includes the core modules, registry queries, and a report section. Thanks for being here from day one!",
    features: [
      "Released app",
      "Launched in VM",
      "Landing page with hero, features, and footer sections",
    ],
  },
  {
    version: "2.0.0",
    date: "2026-08-13",
    title: "Assemblies workflow and transcriptomic registry update",
    description:
      "This release expands the assemblies workflow, adds a dedicated transcriptomic registry page, and improves result handling for filtering, selection, and downloads.",
    features: [
      "Added a dedicated transcriptomic registry page with taxon ID query input and sortable results table",
      "Added deep links from the Assemblies table Transcr. reg. lowest column into the transcriptomic registry page",
      "Added row selection and copy selected GCAs actions to the Assemblies table",
      "Added pipeline filtering to the Assemblies search flow",
        "Added new handover helper for genebuilders",
        "Added longread column to Assemblies table",
    ],
    changes: [
      "Updated the Assemblies results panel with quick filters for short-read and long-read evidence",
      "Added selected GCA copy and full table download actions to Assemblies and transcriptomic registry result tables",
      "Refined table actions and layout around result downloads and selection workflows",
        "Added breed and filtering to the Project pages"
    ],
    fixes: [
      "Added a warning state when Assemblies quick filters hide every row in the table",
      "Reset Assemblies row selection correctly when fetching a new result set",
    ],
  },
 
].reverse();

const formatDate = (date: Date) => {
  return date.toLocaleDateString("en-US", {
    month: "short",
    day: "numeric",
    year: "numeric",
  });
};

export default function Changelog() {
  return (
    <section className="mx-auto max-w-3xl px-6 py-16">
      <h2 className="text-balance font-medium text-4xl tracking-tight">
        Changelog
      </h2>
      <p className="mt-2 text-balance text-lg text-muted-foreground tracking-[-0.015em] sm:mt-3 sm:text-xl">Track all the new features, updates, and fixes in one
        place.
      </p>

      <Separator className="mt-9" />

      <div className="mt-10 flex flex-col gap-6">
        {changelogs.map((changelog) => (
          <div
            className="rounded-xl bg-muted px-6 py-8 sm:px-8"
            key={changelog.version}
          >
            <div>
              <div className="flex items-center text-muted-foreground tracking-tight">
                v{changelog.version} <Dot />{" "}
                {formatDate(new Date(changelog.date))}
              </div>
              <h3 className="mt-3 font-medium text-2xl tracking-[-0.02em]">
                {changelog.title}
              </h3>
              <span className="mt-2 block text-lg text-muted-foreground tracking-tight sm:hidden">
                {changelog.date}
              </span>

              <p className="mt-3 text-foreground/80">{changelog.description}</p>

              <div className="mt-4 space-y-4 text-foreground/80">
                {changelog.features && (
                  <div>
                    <ChangelogSection
                      items={changelog.features}
                      title="Features:"
                    />
                  </div>
                )}

                {changelog.changes && (
                  <ChangelogSection
                    items={changelog.changes}
                    title="Changes:"
                  />
                )}

                {changelog.fixes && (
                  <ChangelogSection items={changelog.fixes} title="Fixes:" />
                )}
              </div>
            </div>
          </div>
        ))}
      </div>
    </section>
  );
}

const ChangelogSection = ({
  title,
  items,
}: {
  title: string;
  items: string[];
}) => {
  return (
    <div>
      <h4 className="mb-1 flex items-center gap-2 font-medium text-foreground text-lg">
        {title}
      </h4>
      <ul className="list-disc pl-5">
        {items.map((item) => (
          <li key={item}>{item}</li>
        ))}
      </ul>
    </div>
  );
};
