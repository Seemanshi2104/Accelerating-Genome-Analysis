class SuffixTreeNode:
    def _init_(self):
        self.children = {}  # key: char, value: (start, length, child node)
        self.suffix_index = -1  # only for leaf nodes

class SuffixTree:
    def _init_(self, text):
        self.text = text
        self.root = SuffixTreeNode()
        self.build_suffix_tree()

    def build_suffix_tree(self):
        n = len(self.text)
        for i in range(n):  # insert all suffixes
            current = self.root
            j = i
            while j < n:
                c = self.text[j]
                if c not in current.children:
                    # create new edge from j to end of string
                    child = SuffixTreeNode()
                    current.children[c] = (j, n - j, child)
                    child.suffix_index = i
                    break
                else:
                    start, length, child = current.children[c]
                    label = self.text[start:start+length]
                    k = 0
                    while k < length and j+k < n and self.text[start+k] == self.text[j+k]:
                        k += 1
                    if k == length:
                        current = child
                        j += k
                    else:
                        # split the edge
                        split = SuffixTreeNode()
                        current.children[c] = (start, k, split)
                        split.children[self.text[start+k]] = (start+k, length-k, child)
                        new_leaf = SuffixTreeNode()
                        split.children[self.text[j+k]] = (j+k, n - (j+k), new_leaf)
                        new_leaf.suffix_index = i
                        break
    def search(self, query):
        current = self.root
        i = 0
        while i < len(query):
            c = query[i]
            if c not in current.children:
                return []  # no match
            start, length, child = current.children[c]
            label = self.text[start:start+length]
            j = 0
            while j < length and i < len(query):
                if query[i] != label[j]:
                    return []  # mismatch
                i += 1
                j += 1
            # even if query ends in middle of an edge, it's a match
            current = child
        # collect all matching suffix indices
        return self.collect_leaf_indices(current)


    def collect_leaf_indices(self, node):
        if node.suffix_index != -1:
            return [node.suffix_index]
        indices = []
        for _, _, child in node.children.values():
            indices.extend(self.collect_leaf_indices(child))
        return indices
    # --- Your Existing Suffix Tree Code (Use exactly as you posted) ---

# Add this after your class definitions:

def assemble_genome_from_reads(reference, reads):
    tree = SuffixTree(reference + "$")
    matched_reads = []

    for read in reads:
        positions = tree.search(read)
        if positions:
            matched_reads.append((positions[0], read))  # use first occurrence
        else:
            print(f"Warning: Read '{read}' not found in reference!")

    # Sort reads based on starting position
    matched_reads.sort()

    # Reconstruct genome using overlap stitching
    assembled = matched_reads[0][1] if matched_reads else ""
    for i in range(1, len(matched_reads)):
        prev = assembled
        _, read = matched_reads[i]

        # Find overlap
        max_overlap = 0
        for j in range(1, min(len(prev), len(read)) + 1):
            if prev[-j:] == read[:j]:
                max_overlap = j
        assembled += read[max_overlap:]  # Add only non-overlapping part

    return assembled

# --- Example Usage ---

reference = "ATCGTTAGCA"
reads = ["TAG", "ATC", "GTT", "AGCA"]

assembled = assemble_genome_from_reads(reference, reads)
print("Assembled genome:", assembled)
