"use client"
import * as React from "react"

import {
    Dialog,
    DialogContent,
    DialogDescription, DialogFooter,
    DialogHeader,
    DialogTitle,
    DialogTrigger,
} from "@/components/ui/dialog"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Button } from "@/components/ui/button"

export function StartAnnotationDialog() {
  return (
    <Dialog>
      <DialogTrigger asChild>
        <Button variant="outline">Start annotation</Button>
      </DialogTrigger>
      <DialogContent className="sm:max-w-[425px]">
        <DialogHeader>
          <DialogTitle>Start annotation</DialogTitle>
          <DialogDescription>
            Generate config file to start annotation.
          </DialogDescription>
        </DialogHeader>
        <div className="grid gap-4 py-4">
          <div className="grid grid-cols-4 items-center gap-4">
            <Label htmlFor="genebuilder_id" className="text-right">
              Genebuilder ID
            </Label>
            <Input id="genebuilder_id" value="60" className="col-span-3" />
          </div>
          <div className="grid grid-cols-4 items-center gap-4">
            <Label htmlFor="workking_dir" className="text-right">
              Working directory
            </Label>
            <Input id="workking_dir" value="working/directory/path" className="col-span-3" />
          </div>
        </div>
        <DialogFooter>
          <Button type="submit">Genertate init file</Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
