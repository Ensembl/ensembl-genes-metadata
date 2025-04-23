import Link from "next/link"

import { MainNav } from "@/components/main-nav"

export function SiteHeader() {
  return (
    <header className="border-grid sticky top-0 z-50 border-b bg-background/95 backdrop-blur supports-[backdrop-filter]:bg-background/60">
      <div className="container-wrap h-20 px-12">
        <div className="flex items-center gap-2 md:gap-4 px-16 py-4">
          <MainNav />
        </div>
      </div>
    </header>
  )
}