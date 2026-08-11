"use client"

import Link from "next/link"
import { Icons } from "@/components/icons"
import { ModeSwitcher } from "@/components/mode_switcher"
import { Search, ListCheck, ChartNoAxesColumnIncreasing, FolderRoot, DatabaseSearch, UserSearch, FlaskConical } from "lucide-react"
import { NavigationMenu,
  NavigationMenuContent,
  NavigationMenuItem,
  NavigationMenuLink,
  NavigationMenuList,
  NavigationMenuTrigger,} from "@/components/ui/navigation-menu"


const sections = [

  {
    title: 'Overview',
    items: [
        {
        icon: (
          <Search />
        ),
        label: 'Find GCA',
        description: 'Search for annotation status',
        href: '/find-gca'
      },
      {
        icon: (
          <ChartNoAxesColumnIncreasing />
        ),
        label: 'Report',
        description: 'Generate reports',
        href: '/report'
      },
      {
        icon: (
          <FolderRoot />
        ),
        label: 'Projects',
        description: 'Quick overview of projects',
        href: '/projects'
      }
    ]
  },
]

const sections_genebuild = [
  {
    title: 'Query the registry',
    items: [
      {
        icon: (
          <DatabaseSearch />
        ),
        label: 'Assemblies',
        description: 'Find out what to annotate',
        href: '/assemblies'
      },
      {
        icon: (
          < ListCheck />
        ),
        label: 'Annotations',
        description: 'See what has been annotated',
        href: '/annotations'
      }
    ]
  },
    {
    title: 'Genebuild info',
    items: [
      {
        icon: (
          <UserSearch />
        ),
        label: 'Dashboard',
        description: 'Personalised dashboard',
        href: '/genebuild'
      },
    ]
  },
 {
    title: 'QC',
    items: [
      {
        icon: (
          <FlaskConical />
        ),
        label: 'Annotation QC',
        description: 'Check annotation quality',
        href: '/annotation-qc'
      },
    ]
  },
]


export function MainNav() {
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


      <NavigationMenu className='*:last:left-auto'>
      <NavigationMenuList>
        <NavigationMenuItem>
          <NavigationMenuTrigger className='[&>svg]:size-4'>Search Genebuild</NavigationMenuTrigger>
          <NavigationMenuContent>
            <div className='w-70 p-2 sm:w-80'>
              <div className='space-y-4'>
                {sections.map(section => (
                  <div key={section.title}>
                    <h4 className='mb-2 text-sm font-medium'>{section.title}</h4>
                    <div className='space-y-1'>
                      {section.items.map(item => (
                        <NavigationMenuLink
                          className='flex items-start gap-2 *:[svg]:mt-1 *:[svg]:size-5'
                          href={item.href}
                          key={item.label}
                        >
                          {item.icon}
                          <div>
                            <p className='font-medium'>{item.label}</p>
                            <p className='text-muted-foreground'>{item.description}</p>
                          </div>
                        </NavigationMenuLink>
                      ))}
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </NavigationMenuContent>
        </NavigationMenuItem>
        <NavigationMenuItem>
          <NavigationMenuTrigger className='[&>svg]:size-4'>Genebuild Tools</NavigationMenuTrigger>
          <NavigationMenuContent>
            <div className='w-70 p-2 sm:w-80'>
              <div className='space-y-4'>
                {sections_genebuild.map(sections_genebuild => (
                  <div key={sections_genebuild.title}>
                    <h4 className='mb-2 text-sm font-medium'>{sections_genebuild.title}</h4>
                    <div className='space-y-1'>
                      {sections_genebuild.items.map(item => (
                        <NavigationMenuLink
                          className='flex items-start gap-2 *:[svg]:mt-1 *:[svg]:size-5'
                          href={item.href}
                          key={item.label}
                        >
                          {item.icon}
                          <div>
                            <p className='font-medium'>{item.label}</p>
                            <p className='text-muted-foreground'>{item.description}</p>
                          </div>
                        </NavigationMenuLink>
                      ))}
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </NavigationMenuContent>
        </NavigationMenuItem>
        <NavigationMenuItem>
        <NavigationMenuLink asChild>
          <Link
            href="/help"
            className="group inline-flex h-10 w-max items-center justify-center rounded-md bg-background px-4 py-2 text-sm font-medium hover:bg-accent hover:text-accent-foreground"
          >
            Help
          </Link>
        </NavigationMenuLink>
      </NavigationMenuItem>
      </NavigationMenuList>
    </NavigationMenu>

      <ModeSwitcher />
    </div>
  </div>
)
}
