"use client"

import * as React from "react"
import { Slot } from "@radix-ui/react-slot"

import { cn } from "@/lib/utils"

// Types
type TimelineContextValue = {
  activeStep: number
  setActiveStep: (step: number) => void
}

// Context
const TimelineContext = React.createContext<TimelineContextValue | undefined>(
  undefined
)

const useTimeline = () => {
  const context = React.useContext(TimelineContext)
  if (!context) {
    throw new Error("useTimeline must be used within a Timeline")
  }
  return context
}

// Timeline
interface TimelineProps extends React.ComponentProps<"div"> {
  defaultValue?: number
  value?: number
  onValueChange?: (value: number) => void
  orientation?: "horizontal" | "vertical"
  asChild?: boolean
}

function Timeline({
  defaultValue = 1,
  value,
  onValueChange,
  orientation = "vertical",
  className,
  asChild,
  ...props
}: TimelineProps) {
  const [activeStep, setInternalStep] = React.useState(defaultValue)

  const setActiveStep = React.useCallback(
    (step: number) => {
      if (value === undefined) {
        setInternalStep(step)
      }
      onValueChange?.(step)
    },
    [value, onValueChange]
  )

  const currentStep = value ?? activeStep
  const Comp = asChild ? Slot : "div"

  return (
    <TimelineContext.Provider
      value={{ activeStep: currentStep, setActiveStep }}
    >
      <Comp
        data-slot="timeline"
        data-orientation={orientation}
        className={cn(
          "group/timeline flex data-[orientation=horizontal]:w-full data-[orientation=horizontal]:flex-row data-[orientation=vertical]:flex-col",
          className
        )}
        {...props}
      />
    </TimelineContext.Provider>
  )
}

// TimelineContent
function TimelineContent({
  className,
  asChild,
  ...props
}: React.ComponentProps<"div"> & { asChild?: boolean }) {
  const Comp = asChild ? Slot : "div"

  return (
    <Comp
      data-slot="timeline-content"
      className={cn("text-muted-foreground text-sm", className)}
      {...props}
    />
  )
}

// TimelineDate
function TimelineDate({
  className,
  asChild,
  ...props
}: React.ComponentProps<"time"> & { asChild?: boolean }) {
  const Comp = asChild ? Slot : "time"

  return (
    <Comp
      data-slot="timeline-date"
      className={cn(
        "mb-1 block font-medium text-muted-foreground text-xs group-data-[orientation=vertical]/timeline:max-sm:h-4",
        className
      )}
      {...props}
    />
  )
}

// TimelineHeader
function TimelineHeader({
  className,
  asChild,
  ...props
}: React.ComponentProps<"div"> & { asChild?: boolean }) {
  const Comp = asChild ? Slot : "div"

  return <Comp data-slot="timeline-header" className={cn(className)} {...props} />
}

// TimelineIndicator
function TimelineIndicator({
  className,
  asChild,
  ...props
}: React.ComponentProps<"div"> & { asChild?: boolean }) {
  const Comp = asChild ? Slot : "div"

  return (
    <Comp
      aria-hidden
      data-slot="timeline-indicator"
      className={cn(
        "group-data-[orientation=horizontal]/timeline:-top-6 group-data-[orientation=horizontal]/timeline:-translate-y-1/2 group-data-[orientation=vertical]/timeline:-left-6 group-data-[orientation=vertical]/timeline:-translate-x-1/2 absolute size-4 rounded-full border-2 border-primary/20 group-data-[orientation=vertical]/timeline:top-0 group-data-[orientation=horizontal]/timeline:left-0 group-data-completed/timeline-item:border-primary",
        className
      )}
      {...props}
    />
  )
}

// TimelineItem
interface TimelineItemProps extends React.ComponentProps<"div"> {
  step: number
  asChild?: boolean
}

function TimelineItem({ step, className, asChild, ...props }: TimelineItemProps) {
  const { activeStep } = useTimeline()
  const Comp = asChild ? Slot : "div"

  return (
    <Comp
      data-slot="timeline-item"
      data-completed={step <= activeStep || undefined}
      className={cn(
        "group/timeline-item relative flex flex-1 flex-col gap-0.5 group-data-[orientation=vertical]/timeline:ms-8 group-data-[orientation=horizontal]/timeline:mt-8 group-data-[orientation=horizontal]/timeline:not-last:pe-8 group-data-[orientation=vertical]/timeline:not-last:pb-6 has-[+[data-completed]]:**:data-[slot=timeline-separator]:bg-primary",
        className
      )}
      {...props}
    />
  )
}

// TimelineSeparator
function TimelineSeparator({
  className,
  asChild,
  ...props
}: React.ComponentProps<"div"> & { asChild?: boolean }) {
  const Comp = asChild ? Slot : "div"

  return (
    <Comp
      aria-hidden
      data-slot="timeline-separator"
      className={cn(
        "group-data-[orientation=horizontal]/timeline:-top-6 group-data-[orientation=horizontal]/timeline:-translate-y-1/2 group-data-[orientation=vertical]/timeline:-left-6 group-data-[orientation=vertical]/timeline:-translate-x-1/2 absolute self-start bg-primary/10 group-last/timeline-item:hidden group-data-[orientation=horizontal]/timeline:h-0.5 group-data-[orientation=vertical]/timeline:h-[calc(100%-1rem-0.25rem)] group-data-[orientation=horizontal]/timeline:w-[calc(100%-1rem-0.25rem)] group-data-[orientation=vertical]/timeline:w-0.5 group-data-[orientation=horizontal]/timeline:translate-x-4.5 group-data-[orientation=vertical]/timeline:translate-y-4.5",
        className
      )}
      {...props}
    />
  )
}

// TimelineTitle
function TimelineTitle({
  className,
  asChild,
  ...props
}: React.ComponentProps<"h3"> & { asChild?: boolean }) {
  const Comp = asChild ? Slot : "h3"

  return (
    <Comp
      data-slot="timeline-title"
      className={cn("font-medium text-sm", className)}
      {...props}
    />
  )
}

export {
  Timeline,
  TimelineContent,
  TimelineDate,
  TimelineHeader,
  TimelineIndicator,
  TimelineItem,
  TimelineSeparator,
  TimelineTitle,
}
