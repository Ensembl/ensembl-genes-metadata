import Link from "next/link"

import { MainNav } from "@/components/main-nav"

export function SiteHeader() {
  return (
    <header className="border-grid sticky top-0 z-50 border-b bg-background/95 backdrop-blur supports-[backdrop-filter]:bg-background/60">
      <div className="min-h-20 w-full max-w-full overflow-visible px-3 sm:px-6 lg:px-12">
        <div className="flex min-w-0 items-center gap-2 py-4 md:gap-4 lg:px-16">
          <MainNav />
        </div>
      </div>
    </header>
  )
}
