"use client";

import React, { useState, useCallback, useRef } from "react";
import { CheckIcon, ChevronsUpDownIcon } from "lucide-react";
import { cn } from "@/lib/utils";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import {
  Command,
  CommandEmpty,
  CommandGroup,
  CommandInput,
  CommandItem,
  CommandList,
} from "@/components/ui/command";
import {
  Popover,
  PopoverContent,
  PopoverTrigger,
} from "@/components/ui/popover";

// Custom debounce hook with proper TypeScript typing
function useDebounceCallback<T extends unknown[]>(
  callback: (...args: T) => void,
  delay: number
) {
  const timeoutRef = useRef<NodeJS.Timeout | null>(null);

  return useCallback(
    (...args: T) => {
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }
      timeoutRef.current = setTimeout(() => {
        callback(...args);
      }, delay);
    },
    [callback, delay]
  );
}

interface TaxonomySearchProps {
  value: string | null;
  onValueChange: (value: string | null) => void;
  className?: string;
}

export default function TaxonomySearch({ value, onValueChange, className }: TaxonomySearchProps) {
  const [open, setOpen] = useState(false);
  const [inputValue, setInputValue] = useState("");
  const [options, setOptions] = useState<string[]>([]);

  const fetchSuggestions = async (query: string): Promise<void> => {
    if (!query) {
      setOptions([]);
      return;
    }
    try {
      const res = await fetch(`http://127.0.0.1:8000/api/taxonomy_search/taxonomy/search?q=${encodeURIComponent(query)}`);
      const data: string[] = await res.json();
      setOptions(data);
    } catch (error) {
      console.error("Failed to fetch suggestions", error);
      setOptions([]);
    }
  };

  // Use custom debounce hook with proper typing
  const debouncedFetch = useDebounceCallback(fetchSuggestions, 300);

  const selectedLabel = value || "Select taxon...";

  return (
    <div className={className}>
      <Label className="mb-2 block">Taxon Search</Label>
      <Popover open={open} onOpenChange={setOpen}>
        <PopoverTrigger asChild>
          <Button
            variant="outline"
            role="combobox"
            aria-expanded={open}
            className="w-[300px] justify-between"
          >
            {selectedLabel}
            <ChevronsUpDownIcon className="ml-2 h-4 w-4 shrink-0 opacity-50" />
          </Button>
        </PopoverTrigger>
        <PopoverContent className="w-[300px] p-0">
          <Command>
            <CommandInput
              placeholder="Type taxon ID or name"
              value={inputValue}
              onValueChange={(val) => {
                setInputValue(val);
                debouncedFetch(val);
              }}
              autoFocus
            />
            <CommandList>
              <CommandEmpty>No taxon found.</CommandEmpty>
              <CommandGroup>
                {Array.isArray(options) &&
                  options.map((option) => (
                    <CommandItem
                      key={option}
                      value={option}
                      onSelect={(currentValue) => {
                        const newValue = currentValue === value ? null : currentValue;
                        onValueChange(newValue);
                        setOpen(false);
                        setInputValue("");
                        setOptions([]);
                      }}
                    >
                      <CheckIcon
                        className={cn(
                          "mr-2 h-4 w-4",
                          value === option ? "opacity-100" : "opacity-0"
                        )}
                      />
                      {option}
                    </CommandItem>
                  ))}
              </CommandGroup>
            </CommandList>
          </Command>
        </PopoverContent>
      </Popover>
    </div>
  );
}