import { SiGithub } from "@icons-pack/react-simple-icons";
import Link from "next/link";
import {Icons} from "@/components/icons";

export default function Footer() {
  return (
    <footer className="border-t bg-background px-6 py-3">
      <div className="mx-auto flex w-full max-w-screen-2xl items-center justify-between">
        <Link href="/" className="flex items-center gap-2">
          <Icons.logo className="h-5 w-5" />
          <span className="text-sm">Genebuild Metadata</span>
        </Link>

        <div className="flex items-center gap-6 text-sm">
          <Link
            href="/help"
            target="_blank"
            rel="noopener noreferrer"
            className="hover:underline"
          >
            Help
          </Link>
          <Link
            href="http://genebuild-metadata.ebi.ac.uk:8000/docs"
            target="_blank"
            rel="noopener noreferrer"
            className="hover:underline"
          >
            API Docs
          </Link>

          <Link
            href="https://github.com/Ensembl/ensembl-genes-metadata/tree/dev/gb_metadata_handling/metadata_app"
            target="_blank"
            rel="noopener noreferrer"
            className="text-muted-foreground hover:text-foreground"
          >
            <SiGithub className="h-5 w-5" />
          </Link>
        </div>
      </div>
    </footer>
  );
}
