"use client";
import React, { useState } from "react"
import { CardsStats } from "@/components/ui/cards_stats"
import { CardsCalendar } from "@/components/ui/card_calendar"
import  { CardsDataTable } from "@/components/ui/card_projects"
import { WelcomeCard } from "@/components/ui/card_welcome"

export default function Page() {
  return (
      <div className="min-h-screen justify-center flex flex-wrap align-items-center">
          <div className="container m-16 mt-10 max-w-7xl">
            <div className="pt-8 pb-10 py-16">
              <div className="grid grid-cols-1 md:grid-cols-1 gap-6 lg:grid-cols-2 justify-items-center">
                  <div className="w-full">
                    <WelcomeCard></WelcomeCard>
                  </div>
                  <div className="w-full">
                    <CardsStats></CardsStats>
                  </div>
                  <div className="w-full">
                    <CardsCalendar></CardsCalendar>
                  </div>
                 <div className="w-full">
                    <CardsDataTable></CardsDataTable>
                 </div>
              </div>
            </div>
          </div>
      </div>
  )
}