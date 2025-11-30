def test_extension_length_limit(self):
    blast_101_search.word_size = 3
    blast_101_search.tdict.word_size = 3

    myline_database = "ATGC"
    blast_101_search.args.query = "ATGC"
    blast_101_search.qsequence = blast_101_search.args.query
    blast_101_search.query_sequence = blast_101_search.tdict.create_word_dict(
        blast_101_search.qsequence
    )

    # Patch programme_settings (if needed)
    with patch.dict(blast_101_search.programme_settings.settings["BLAST"], {
        "min_extension_score": "0",
        "max_extension_length": "10",
    }):
        bestscore = blast_101_search.process_blast(myline_database)
        print(bestscore)
