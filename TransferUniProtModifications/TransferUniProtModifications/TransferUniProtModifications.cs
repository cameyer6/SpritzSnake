using Fclp;
using System;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TransferUniProtModifications
{
    internal class TransferUniProtModifications
    {
        private static void Main(string[] args)
        {
            Console.WriteLine("Welcome to TransferModifications!");
            var p = new FluentCommandLineParser<ApplicationArguments>();

            p.Setup(arg => arg.UniProtXml)
                .As('x', "uniprot_xml")
                .Required()
                .WithDescription("UniProt protein XML file.");

            p.Setup(arg => arg.SpritzXml)
                .As('y', "spritz_xml")
                .Required()
                .WithDescription("Custom protein XML file, e.g. from Spritz.");

            p.SetupHelp("h", "help")
                .Callback(text => Console.WriteLine(text));

            var result = p.Parse(args);

            TransferModifications(p.Object.UniProtXml, p.Object.SpritzXml);
        }

        public static string TransferModifications(string sourceXmlPath, string destinationXmlPath)
        {
            var uniprotPtms = ProteinAnnotation.GetUniProtMods(Environment.CurrentDirectory);
            var uniprot = ProteinDbLoader.LoadProteinXML(sourceXmlPath, true, DecoyType.None, uniprotPtms, false, null, out var un);
            string outxml = Path.Combine(Path.GetDirectoryName(destinationXmlPath), Path.GetFileNameWithoutExtension(destinationXmlPath) + ".withmods.xml");
            var nonVariantProts = ProteinDbLoader.LoadProteinXML(destinationXmlPath, true, DecoyType.None, uniprotPtms, false, null, out un).Select(p => p.NonVariantProtein).Distinct();
            var newProts = ProteinAnnotation.CombineAndAnnotateProteins(uniprot, nonVariantProts.ToList());
            ProteinDbWriter.WriteXmlDatabase(null, newProts, outxml);
            string outfasta = Path.Combine(Path.GetDirectoryName(destinationXmlPath), Path.GetFileNameWithoutExtension(destinationXmlPath) + ".fasta");
            var prot = newProts.FirstOrDefault(p => p.Accession.Contains("_"));
            ProteinDbWriter.WriteFastaDatabase(newProts.SelectMany(p => p.GetVariantProteins()).ToList(), outfasta, "|");
            return outxml;
        }

        public class ApplicationArguments
        {
            public string SpritzXml { get; set; }
            public string ReferenceGeneModel { get; set; }
            public string UniProtXml { get; set; }
        }
    }
}