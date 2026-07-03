"use client"

import Link from "next/link"
import { usePathname } from "next/navigation"
import { cn } from "@/lib/utils"
import { Icons } from "@/components/icons"
import { ModeSwitcher } from "@/components/mode_switcher"

export function MainNav() {
  const pathname = usePathname()

  return (
  <div className="flex w-full min-w-0 flex-wrap items-center gap-y-3 sm:flex-nowrap">
    {/* Logo Section */}
    <Link href="/" className="mr-4 flex shrink-0 items-start gap-2 lg:mr-6">
      <Icons.logo className="h-6 w-6" />
      <span className="hidden font-bold lg:inline-block leading-tight">
        Genebuild<br />Metadata
      </span>
    </Link>

    {/* Right-side nav + mode switch */}
    <div className="ml-auto flex min-w-0 flex-1 flex-wrap items-center justify-end gap-3 sm:flex-nowrap sm:gap-6">
      <nav className="flex min-w-0 flex-wrap items-center justify-end gap-x-3 gap-y-2 overflow-visible text-sm sm:flex-nowrap sm:gap-4 xl:gap-6">
        <Link
          href="/"
          className={cn(
            "transition-colors hover:text-foreground/80",
            pathname === "/"
              ? "text-foreground font-bold"
              : "text-foreground/80"
          )}
        >
          Home
        </Link>
        <Link
          href="/assemblies"
          className={cn(
            "transition-colors hover:text-foreground/80",
            pathname?.startsWith("/assemblies")
              ? "text-foreground font-bold"
              : "text-foreground/80"
          )}
        >
          Assemblies
        </Link>
        <Link
          href="/annotations"
          className={cn(
            "transition-colors hover:text-foreground/80",
            pathname?.startsWith("/annotations")
              ? "text-foreground font-bold"
              : "text-foreground/80"
          )}
        >
          Annotations
        </Link>
        <Link
          href="/report"
          className={cn(
            "transition-colors hover:text-foreground/80",
            pathname?.startsWith("/report")
              ? "text-foreground font-bold"
              : "text-foreground/80"
          )}
        >
          Report
        </Link>
        <Link
          href="/projects"
          className={cn(
            "transition-colors hover:text-foreground/80",
            pathname?.startsWith("/projects")
              ? "text-foreground font-bold"
              : "text-foreground/80"
          )}
        >
          Projects
        </Link>
        <Link
          href="/genebuild"
          className={cn(
            "transition-colors hover:text-foreground/80",
            pathname?.startsWith("/genebuild")
              ? "text-foreground font-bold"
              : "text-foreground/80"
          )}
        >
          Genebuild
        </Link>
      </nav>

      <ModeSwitcher />
    </div>
  </div>
)
}
