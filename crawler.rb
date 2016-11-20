require 'rubygems'
require 'nokogiri'
require 'open-uri'

DATA_DIR = "fastas"
Dir.mkdir(DATA_DIR) unless File.exists?(DATA_DIR)

base_url = "http://pax-db.org/protein/"

# load ids and abundances
molecules = []
infile = File.open("9606-PA_plasma.txt","r")
while line = infile.gets
	unless line[0] == '#' # header
		parts = line.chop.split(" ")
		molecules << [parts[0],parts[2]]
	end
end
infile.close

outfile = File.open("9606-PA_plasma.fasta","w")
log = File.open("9606-PA_plasma.log","w")

molecules.each_w_index do |molecule, m_idx|
	puts "Building #{m_idx} of #{molecules.size}"

	#go to pax-db page for this molecule
	page = Nokogiri::HTML(open("#{base_url}#{molecule[0]}"))

	#get uniprot page link
	link = page.xpath("/html/body/div[2]/p/a")[0].attributes["href"].value

	#pull fasta entry from uniprot
	uniprot_molecule_page = ""
	begin
	uniprot_molecule_page = Nokogiri::HTML(open(link))

	fasta_link = uniprot_molecule_page.css("//div#sequences.section div#sequences-section div.sequence-isoform div.sequence-isoform-leftcol a.tooltipped.icon.icon-functional.button.displayThis")[0].attributes["href"]

	fasta = Nokogiri::HTML(open("http://uniprot.org#{fasta_link}")).elements[0].child.child.child.text

	#write out fasta file, inserting standard abundance
	fasta_lines = fasta.split("\n")
	outfile.puts "#{fasta_lines[0]} ##{molecule[1]}" 
	fasta_lines[1..-1].each do |fasta_line|
		outfile.puts fasta_line
	end
	rescue
		log.puts "Failed to open #{link}"
	end
end

outfile.close
log.close
