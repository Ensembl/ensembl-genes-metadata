"use client"

import Link from "next/link"
import { usePathname } from "next/navigation"
import { cn } from "@/lib/utils"
import { Icons } from "@/components/icons"
import { ModeSwitcher } from "@/components/mode_switcher"

export function MainNav() {
  const pathname = usePathname()

  return (
  <div className="flex w-full items-center">
    {/* Logo Section */}
    <Link href="/" className="flex items-start mr-6 gap-2 lg:mr-6">
      <Icons.logo className="h-6 w-6" />
      <span className="hidden font-bold lg:inline-block leading-tight">
        Genebuild<br />Metadata
      </span>
    </Link>

    {/* Right-side nav + mode switch */}
    <div className="flex items-center gap-6 ml-auto">
      <nav className="flex items-center gap-4 text-sm xl:gap-6">
        <Link
          href="/"
          className={cn(
            "transition-colors hover:text-foreground/80",
            pathname === "/"
              ? "text-foreground"
              : "text-foreground/80"
          )}
        >
          Home
        </Link>
        <Link
          href="/assemblies"
          className={cn(
            "transition-colors hover:text-foreground/80",
            pathname === "/assemblies"
              ? "text-foreground"
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
              ? "text-foreground"
              : "text-foreground/80"
          )}
        >
          Annotations
        </Link>
        <Link
          href="/transcriptomic"
          className={cn(
            "transition-colors hover:text-foreground/80",
            pathname?.startsWith("/transcriptomic")
              ? "text-foreground"
              : "text-foreground/80"
          )}
        >
          Transcriptomic assessment
        </Link>
      </nav>

      <ModeSwitcher />
    </div>
  </div>
)
}