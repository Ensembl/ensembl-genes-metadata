"use client"

import * as React from "react"
import { Button } from "@/components/ui/button"
import { Popover, PopoverTrigger, PopoverContent } from "@/components/ui/popover"
import { CheckIcon, SlidersHorizontal, ListFilterPlus } from "lucide-react"
import { cn } from "@/lib/utils"

function PopoverCheckboxItem({
  checked,
  onCheckedChange,
  children,
  disabled = false,
  ...props
}: {
  checked: boolean
  onCheckedChange: () => void
  children: React.ReactNode
  disabled?: boolean
}) {
  return (
    <div
      className={cn(
        "flex items-center gap-2 select-none p-2 rounded-sm text-sm",
        disabled ? "text-muted-foreground cursor-not-allowed" : "hover:bg-accent"
      )}
      onClick={disabled ? undefined : onCheckedChange}
      {...props}
    >
      <div
        className={cn(
          "w-4 h-4 flex justify-center items-center",
          checked ? "text-primary" : "bg-transparent"
        )}
      >
        {checked && <CheckIcon className="size-4" />}
      </div>
      <span>{children}</span>
    </div>
  )
}

interface PopoverWithMultiSelectProps {
  selectedItems: string[]
  setSelectedItems: React.Dispatch<React.SetStateAction<string[]>>
  onAutoFillHighQuality: () => void
}

export function PopoverWithMultiSelect({
  selectedItems,
  setSelectedItems,
  onAutoFillHighQuality,
}: PopoverWithMultiSelectProps)  {

  const handleToggleItem = (item: string) => {
    setSelectedItems((prevSelectedItems) =>
      prevSelectedItems.includes(item)
        ? prevSelectedItems.filter((i) => i !== item)
        : [...prevSelectedItems, item]
    )
  }


  const group1 = ["Assembly level", "Assembly type", "Contig N50", "Sequence length"];
  const group2 = ["GC%", "Genome coverage", "Number of contigs", "Number of scaffolds", "Scaffold N50"];

  return (
    <Popover>
      <PopoverTrigger asChild>
        <Button variant="outline" className="w-full dark:bg-foreground dark:text-background">
          <ListFilterPlus className="shrink-0" />
          <span className="hidden lg:inline">Assembly Metrics</span>
        </Button>
      </PopoverTrigger>

      <PopoverContent className="w-full">
        <Button variant="ghost" className="w-full justify-start" onClick={onAutoFillHighQuality}>
          <SlidersHorizontal className="" />
          Get annotation candidates
        </Button>

        <div className="bg-border my-2 -mx-2 h-px w-auto" />

        <div className="space-y-2">
          <div className="grid gap-1">
            {group1.map((item) => (
              <PopoverCheckboxItem
                key={item}
                checked={selectedItems.includes(item)}
                onCheckedChange={() => handleToggleItem(item)}
              >
                {item}
              </PopoverCheckboxItem>
            ))}
          </div>

          <div className="bg-border my-2 -mx-2 h-px w-auto" />

          <div className="grid gap-1">
            {group2.map((item) => (
              <PopoverCheckboxItem
                key={item}
                checked={selectedItems.includes(item)}
                onCheckedChange={() => handleToggleItem(item)}
              >
                {item}
              </PopoverCheckboxItem>
            ))}
          </div>
        </div>
      </PopoverContent>
    </Popover>
  );
}