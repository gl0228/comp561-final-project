# Pseudocode

# Run BLAST with specified parameters
blast_results = run_blastn(
    query=query_path,
    database=blast_db_path,
    evalue=1000,
    word_size=word_size,
    perc_identity=0,
    outfmt=6,
    ungapped=True,
    strand="plus"
)

match_probability = empty_map()
max_prob = -∞
max_match = (0, 0)
epsilon = 1e-10

for each row in blast_results:
    probability = log(1)  # This is effectively 0.
    unlikely_match = false

    # Extract indices from BLAST output (convert to zero-based)
    query_start = row.query_start - 1
    query_end = row.query_end - 1
    db_start = row.db_start - 1
    db_end = row.db_end - 1

    diff = db_start - query_start

    # Determine how far to extend to match full query length
    left_extension = min(query_start, db_start)
    right_extension = min(query_length - query_end, database_length - db_end)

    extended_db_start = db_start - left_extension
    extended_db_end = db_end + right_extension - 1
    extended_query_start = query_start - left_extension
    extended_query_end = query_end + right_extension - 1

    start_end_tuple = (extended_db_start, extended_db_end)
    consecutive_low_prob_count = 0

    # Compute probability over the extended range
    for index from extended_db_start to extended_db_end:
        nucleotide = sequence_array[extended_query_start]
        nucleotide_prob = positional_probabilities[index][nucleotide]

        if heuristic is enabled:
            if nucleotide_prob < heuristic_prob:
                consecutive_low_prob_count += 1
            else:
                consecutive_low_prob_count = 0

            if consecutive_low_prob_count >= heuristic_num:
                # This match is considered unlikely; skip it
                unlikely_match = true
                break

        # Add log probability of the nucleotide match
        if nucleotide_prob > 0:
            probability += log(nucleotide_prob)
        else:
            probability += log(epsilon)

        extended_query_start += 1

    if unlikely_match:
        # Do not record this match
        continue

    # Update max probability and store result
    if probability > max_prob:
        max_prob = probability
        max_match = start_end_tuple

    match_probability[start_end_tuple] = probability

# Sort matches by probability in descending order
sorted_match_probabilities = sort_map_by_value(match_probability, descending=true)

return sorted_match_probabilities
