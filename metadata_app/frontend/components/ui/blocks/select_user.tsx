import { Avatar, AvatarFallback, AvatarImage } from "@/components/ui/avatar";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuLabel,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu";
import { Check, ChevronsUpDown } from "lucide-react";
import { useState } from "react";

const people = [
  {
    name: "Vianey",
    role: "vianey",
      imageUrl: "https://content.embl.org/sites/default/files/persons/CP-60041730.jpg",

  },
  {
    name: "Francesca",
    role: "ftricomi",
      imageUrl: "https://content.embl.org/sites/default/files/persons/CP-60035116.jpg",

  },
  {
    name: "Swati",
    role: "swati",
      imageUrl: "https://content.embl.org/sites/default/files/persons/CP-60039720.jpg",

  },
  {
    name: "Jose",
    role: "ereboperezsilva",
    imageUrl: "https://content.embl.org/sites/default/files/persons/CP-60032924.jpg",
  },
      {
    name: "Jack",
    role: "jackt",
          imageUrl: "https://content.embl.org/sites/default/files/persons/CP-60041119.jpg",

  },
      {
    name: "Anna",
    role: "lazar",
          imageUrl: "https://content.embl.org/sites/default/files/persons/CP-60045345.jpg",

  },
    {
    name: "Leanne",
    role: "leanne",
        imageUrl: "https://content.embl.org/sites/default/files/persons/CP-60025510.jpg",

  },
];

type SelectUserProps = {
  onSelect: (user: {
    name: string;
    role: string;
    imageUrl: string;
  }) => void;
  placeholder?: string;
};

export default function GridList02({ onSelect }: SelectUserProps) {
      const [selectedUser, setSelectedUser] = useState<null | typeof people[0]>(null);


  return (
    <DropdownMenu>
      <DropdownMenuTrigger className="flex items-center gap-2 bg-accent py-2.5 px-3 rounded-xl">
        <Avatar className="rounded-xl h-8 w-8">
          {selectedUser ? (
            <AvatarImage src={selectedUser.imageUrl} alt={selectedUser.name} />
          ) : (
            <AvatarFallback className="rounded-xl bg-primary/20 text-primary-foreground">
              ?
            </AvatarFallback>
          )}
        </Avatar>

        <div className="text-start flex flex-col gap-1 leading-none">
          <span className="text-sm font-semibold truncate max-w-[17ch]">
            {selectedUser ? selectedUser.name : "Select a Genebuilder"}
          </span>
          <span className="text-xs text-muted-foreground truncate max-w-[20ch]">
            {selectedUser ? selectedUser.role : ""}
          </span>
        </div>

        <ChevronsUpDown className="ml-6 h-4 w-4 text-muted-foreground" />
      </DropdownMenuTrigger>

      <DropdownMenuContent className="w-52" align="start">
        <DropdownMenuLabel>Genebuilders</DropdownMenuLabel>
        {people.map((person) => (
          <DropdownMenuItem
            key={person.name}
            onClick={() => {
              setSelectedUser(person);
              onSelect(person);
            }}
          >
            <div className="flex items-center gap-2">
              <Avatar className="rounded-md h-8 w-8">
                <AvatarImage src={person.imageUrl} alt={person.name} />
                <AvatarFallback className="rounded-md bg-primary/10 text-foreground">
                  {person.name[0]}
                </AvatarFallback>
              </Avatar>
              <div className="flex flex-col">
                <span>{person.name}</span>
                <span className="text-xs text-muted-foreground">{person.role}</span>
              </div>
            </div>
            {selectedUser?.name === person.name && <Check className="ml-auto" />}
          </DropdownMenuItem>
        ))}
      </DropdownMenuContent>
    </DropdownMenu>
  );
}

