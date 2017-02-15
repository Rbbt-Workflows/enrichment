require 'rbbt'
require 'rbbt/util/misc'
require 'rbbt/workflow'
require 'rbbt/association'
require 'rbbt/knowledge_base'
require 'rbbt/statistics/hypergeometric'
require 'rbbt/statistics/random_walk'

Workflow.require_workflow 'Genomics'
require 'rbbt/knowledge_base/Genomics'

Workflow.require_workflow 'Translation'

module Enrichment
  extend Workflow
  extend Resource

  self.subdir = "var/Enrichment"

  class << self
    attr_accessor :knowledge_base_dir

    def knowledge_base_dir
      @knowledge_base_dir ||= Enrichment.knowledge_base
    end
  end


  MASKED_TERMS = %w(cancer melanoma carcinoma glioma hepatitis leukemia leukaemia disease infection opathy hepatitis sclerosis hepatatis glioma Shigellosis)
  MASKED_IDS = {}

  RENAMES = Organism.identifiers(Organism.default_code("Hsa")).tsv(:persist => false, :key_field => "Ensembl Gene ID", :fields => [], :grep => '^#\|PCDH', :type => :list, :persit_update => true).add_field("Cluster"){ "Cadherin" }
  RENAMES.type = :single
  RENAMES.process("Cluster") do |values|
    values.first
  end

  begin
    RENAMES.keys.to("KEGG Gene ID").compact.each do |gene|
      RENAMES[gene] = "Cadherin"
    end
  rescue
    Log.warn "Could not process KEGG renames"
  end

  Organism.identifiers(Organism.default_code("Hsa")).tsv(:persist => false, :key_field => "UniProt/SwissProt Accession", :fields => [], :grep => '^#\|PCDH', :type => :list, :persit_update => true).add_field("Cluster"){ "Cadherin" }.keys.each do |gene|
    RENAMES[gene] = "Cadherin"
  end

  Organism.identifiers(Organism.default_code("Hsa")).tsv(:persist => false, :key_field => "Entrez Gene ID", :fields => [], :grep => '^#\|PCDH', :type => :list, :persit_update => true).add_field("Cluster"){ "Cadherin" }.keys.each do |gene|
    RENAMES[gene] = "Cadherin"
  end

  REGISTRY = Genomics.knowledge_base.registry
  DATABASES = REGISTRY.keys

  helper :database_info do |database, organism|
    Persist.memory([database, organism] * ": ") do
      @organism_kb ||= {}
      @organism_kb[organism] ||= begin
                                  dir = Enrichment.knowledge_base_dir

                                  kb = KnowledgeBase.new dir, organism
                                  kb.format["Gene"] = "Ensembl Gene ID"
                                  kb.registry = REGISTRY
                                  DATABASES.replace(REGISTRY.keys)
                                  kb
                                end

      db = @organism_kb[organism].get_database(database, :persist => true, :target => "Gene=>Ensembl Gene ID" )
      db = Association.add_reciprocal db if @organism_kb[organism].registry[database][1][:undirected]

      tsv, total_keys, source_field, target_field = [db, db.keys, db.key_field, db.fields.first]

      tsv.filename = db.persistence_path

      if target_field == "Ensembl Gene ID"
        pathway_field, gene_field = source_field, target_field
        total_genes = Gene.setup(tsv.values.flatten.compact.uniq, "Ensembl Gene ID", organism)
      else
        pathway_field, gene_field = target_field, source_field
        total_genes = total_keys
      end

      tsv.namespace = organism

      [tsv, total_genes, gene_field, pathway_field]
    end
  end

  input :database, :select, "Database code", nil, :select_options => DATABASES
  input :list, :array, "Gene list in any supported format; they will be translated accordingly"
  input :organism, :string, "Organism code (not used for kegg)", Organism.default_code("Hsa")
  input :cutoff, :float, "Cutoff value", 0.05
  input :fdr, :boolean, "Perform Benjamini-Hochberg FDR correction", true
  input :background, :array, "Enrichment background", nil
  input :invert_background, :boolean, "Restrict to elements NOT in background", false
  input :mask_diseases, :boolean, "Mask disease related terms", true
  input :fix_clusters, :boolean, "Fixed dependence in gene clusters", true
  task :enrichment => :tsv do |database, list, organism, cutoff, fdr, background, invert_background, mask_diseases, fix_clusters|
    raise ParameterException, "No list given" if list.nil?

    background = nil if Array === background and background.empty?
    ensembl    = Translation.job(:translate, nil, :format => "Ensembl Gene ID", :genes => list, :organism => organism).run.compact.uniq
    background = Translation.job(:translate, nil, :format => "Ensembl Gene ID", :genes => background, :organism => organism).run.compact.uniq if background and background.any?

    ensembl.reject!{|e| e.nil? or e.empty?}

    Gene.setup(ensembl, "Ensembl Gene ID", Organism.default_code("Hsa"))
    Gene.setup(background, "Ensembl Gene ID", Organism.default_code("Hsa")) if background

    database_tsv, all_db_genes, database_key_field, database_field = database_info database, organism

    if invert_background and background
      background = all_db_genes - background
    end

    if background.nil?
      background = all_db_genes - Organism.blacklist_genes(organism).list
    end

    database_tsv.with_unnamed do
      log :reordering, "Reordering database"
      database_tsv.with_monitor :desc => "Reordering" do
        database_tsv = database_tsv.reorder "Ensembl Gene ID", nil, :persist => true, :zipped => true, :merge => true
      end
    end unless "Ensembl Gene ID" == database_field

    if mask_diseases and not Gene == Entity.formats[database_field]
      Log.debug("Masking #{MASKED_TERMS * ", "}")
      masked = MASKED_IDS[database] ||= database_tsv.with_unnamed do
        terms = database_tsv.values.flatten.uniq
        terms = Misc.prepare_entity(terms, database_field)
        if terms.respond_to? :name
          terms.select{|t| t.name =~ /#{MASKED_TERMS * "|"}/i}
        else
          masked = nil
        end
      end
    else
      masked = nil
    end

    database_tsv = database_tsv.to_flat unless database_tsv.type == :flat

    log :enrichment, "Calculating Enrichment"
    database_tsv.enrichment(ensembl, database_field, :persist => (background.nil? or background.empty?), :cutoff => cutoff, :fdr => fdr, :background => background, :rename => (fix_clusters ? Enrichment::RENAMES : nil), :masked => masked).tap{|tsv| tsv.namespace = organism}
  end
  export_synchronous :enrichment

  input :database, :select, "Database code", nil, :select_options => DATABASES
  input :list, :array, "Gene list in any supported format; they will be translated accordingly"
  input :organism, :string, "Organism code (not used for kegg)", Organism.default_code("Hsa")
  input :permutations, :integer, "Number of permutations used to compute p.value", 10000
  input :cutoff, :float, "Cutoff value", 0.05
  input :fdr, :boolean, "Perform Benjamini-Hochberg FDR correction", true
  input :background, :array, "Enrichment background", nil
  input :mask_diseases, :boolean, "Mask disease related terms", true
  input :fix_clusters, :boolean, "Fixed dependence in gene clusters", true
  input :count_missing, :boolean, "Account for genes with pathway annotations that are missing in list", false
  task :rank_enrichment => :tsv do |database, list, organism, permutations, cutoff, fdr, background, mask_diseases, fix_clusters, count_missing|
    raise ParameterException, "No list given" if list.nil?

    ensembl    = Translation.job(:translate, nil, :format => "Ensembl Gene ID", :genes => list, :organism => organism).run.compact.uniq
    background = Translation.job(:translate, nil, :format => "Ensembl Gene ID", :genes => background, :organism => organism).run.compact.uniq if background and background.any?

    ensembl.reject!{|e| e.nil? or e.empty?}

    Gene.setup(ensembl, "Ensembl Gene ID", Organism.default_code("Hsa"))
    Gene.setup(background, "Ensembl Gene ID", Organism.default_code("Hsa")) if background

    database_tsv, all_db_genes, database_key_field, database_field = database_info database, organism

    if mask_diseases and not Gene == Entity.formats[database_field]
      Log.debug("Masking #{MASKED_TERMS * ", "}")
      masked = MASKED_IDS[database] ||= database_tsv.with_unnamed do
        terms = database_tsv.values.flatten.uniq
        terms = Misc.prepare_entity(terms, database_field)
        if terms.respond_to? :name
          terms.select{|t| t.name =~ /#{MASKED_TERMS * "|"}/i}
        else
          masked = nil
        end
      end
    else
      masked = nil
    end

    missing = (all_db_genes - list).length if count_missing

    cutoff = cutoff.to_f

    log :enrichment, "Performing enrichment"
    database_tsv = database_tsv.to_flat unless database_tsv.type == :flat
    database_tsv.rank_enrichment(ensembl,  :persist => (background.nil? or background.empty?), :cutoff => cutoff, :fdr => fdr, :background => background, :rename => (fix_clusters ? RENAMES : nil), :permutations => permutations, :persist_permutations => true, :missing => missing || 0, :masked => masked).select("p-value"){|p| p.to_f.abs <= cutoff}.tap{|tsv| tsv.namespace = organism}
  end
  export_asynchronous :rank_enrichment
end
