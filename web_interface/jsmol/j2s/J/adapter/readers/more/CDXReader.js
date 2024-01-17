Clazz.declarePackage ("J.adapter.readers.more");
Clazz.load (["J.adapter.readers.xml.XmlChemDrawReader"], "J.adapter.readers.more.CDXReader", ["JU.Rdr", "J.adapter.writers.CDXMLWriter"], function () {
c$ = Clazz.declareType (J.adapter.readers.more, "CDXReader", J.adapter.readers.xml.XmlChemDrawReader);
Clazz.defineMethod (c$, "processXml2", 
function (parent, saxReader) {
var xml = J.adapter.writers.CDXMLWriter.fromCDX (this.binaryDoc);
this.reader = JU.Rdr.getBR (xml);
Clazz.superCall (this, J.adapter.readers.more.CDXReader, "processXml2", [this, saxReader]);
this.binaryDoc = null;
}, "J.adapter.readers.xml.XmlReader,~O");
});
