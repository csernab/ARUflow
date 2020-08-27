rule clean:
    shell:
        """
        rm -rf output_data/*
        rm -rf results/*
        rm -r .snakemake/
	rm -r logs/*
        """
