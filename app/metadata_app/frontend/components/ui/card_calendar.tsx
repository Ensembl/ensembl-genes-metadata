"use client"

import { addDays } from "date-fns"

import { Calendar } from "@/components/ui/calendar"
import {Card, CardHeader, CardContent, CardTitle} from "@/components/ui/card"

const start = new Date(2023, 5, 5)

export function CardsCalendar() {
  return (
    <div className="grid grid-cols-2 gap-4">
      <Card className="justify-items-center">
        <CardHeader>
            <CardTitle className="text-lg font-bold">Metadata registry updates</CardTitle>
        </CardHeader>
        <CardContent className="justify-items-center w-full">
            <Calendar
              numberOfMonths={1}
              mode="range"
              defaultMonth={start}
              selected={{
                from: start,
                to: addDays(start, 8),
              }}
            />
          </CardContent>
        </Card>
        <Card className="justify-items-center">
        <CardHeader>
            <CardTitle className="text-lg font-bold">Transcriptomic registry updates</CardTitle>
        </CardHeader>
        <CardContent className="justify-items-center w-full">
            <Calendar
              numberOfMonths={1}
              mode="range"
              defaultMonth={start}
              selected={{
                from: start,
                to: addDays(start, 8),
              }}
            />
          </CardContent>
        </Card>
    </div>
  )
}