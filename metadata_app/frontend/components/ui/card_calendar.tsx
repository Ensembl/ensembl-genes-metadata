"use client"

import { useEffect, useState } from "react"
import { parseISO } from "date-fns"

import { Calendar } from "@/components/ui/calendar"
import { Card, CardHeader, CardContent, CardTitle } from "@/components/ui/card"

async function fetchUpdateDates(endpoint: string): Promise<Date[]> {
  const res = await fetch(endpoint)
  const data = await res.json()
  return data.map((d: string) => parseISO(d)) // Parse ISO strings to Date objects
}

export function CardsCalendar() {
  const [metadataDates, setMetadataDates] = useState<Date[]>([])
  const [transcriptomicDates, setTranscriptomicDates] = useState<Date[]>([])

  useEffect(() => {
    fetchUpdateDates("/api/home_page/home/meta_update").then(setMetadataDates)
    fetchUpdateDates("/api/home_page/home/transc_update").then(setTranscriptomicDates)
  }, [])

  const getLastDate = (dates: Date[]) => {
    if (dates.length > 0) {
      return dates.reduce((latest, current) => (current > latest ? current : latest))
    }
    return new Date()
  }

  return (
    <div className="grid grid-cols-2 gap-4">
      <Card className="dark:bg-secondary">
        <CardHeader>
          <CardTitle className="text-lg font-bold">Metadata registry updates</CardTitle>
        </CardHeader>
        <CardContent className="w-full justify-items-center">
          {metadataDates.length > 0 && (
            <Calendar
              numberOfMonths={1}
              mode="multiple"
              selected={metadataDates}
              defaultMonth={getLastDate(metadataDates)}
            />
          )}
        </CardContent>
      </Card>

      <Card className="dark:bg-secondary">
        <CardHeader>
          <CardTitle className="text-lg font-bold">Transcriptomic registry updates</CardTitle>
        </CardHeader>
        <CardContent className="w-full justify-items-center">
          {transcriptomicDates.length > 0 && (
            <Calendar
              numberOfMonths={1}
              mode="multiple"
              selected={transcriptomicDates}
              defaultMonth={getLastDate(transcriptomicDates)}
            />
          )}
        </CardContent>
      </Card>
    </div>
  )
}