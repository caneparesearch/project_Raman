Clazz.declarePackage ("J.adapter.readers.xml");
Clazz.load (["J.adapter.readers.xml.XmlReader", "J.adapter.smarter.Atom", "$.Bond", "java.util.ArrayList", "$.Stack"], "J.adapter.readers.xml.XmlChemDrawReader", ["java.lang.Boolean", "JU.BS", "$.Lst", "$.PT", "J.api.JmolAdapter", "JU.Edge", "$.Logger"], function () {
c$ = Clazz.decorateAsClass (function () {
this.minX = 3.4028235E38;
this.minY = 3.4028235E38;
this.minZ = 3.4028235E38;
this.maxZ = -3.4028235E38;
this.maxY = -3.4028235E38;
this.maxX = -3.4028235E38;
this.no3D = false;
if (!Clazz.isClassDefined ("J.adapter.readers.xml.XmlChemDrawReader.CDNode")) {
J.adapter.readers.xml.XmlChemDrawReader.$XmlChemDrawReader$CDNode$ ();
}
if (!Clazz.isClassDefined ("J.adapter.readers.xml.XmlChemDrawReader.CDBond")) {
J.adapter.readers.xml.XmlChemDrawReader.$XmlChemDrawReader$CDBond$ ();
}
this.fragments = null;
this.thisFragment = null;
this.thisNode = null;
this.nodes = null;
this.nostereo = null;
this.textBuffer = null;
Clazz.instantialize (this, arguments);
}, J.adapter.readers.xml, "XmlChemDrawReader", J.adapter.readers.xml.XmlReader);
Clazz.prepareFields (c$, function () {
this.fragments =  new java.util.Stack ();
this.nodes =  new java.util.Stack ();
this.nostereo =  new java.util.ArrayList ();
});
Clazz.makeConstructor (c$, 
function () {
Clazz.superConstructor (this, J.adapter.readers.xml.XmlChemDrawReader, []);
});
Clazz.overrideMethod (c$, "processXml", 
function (parent, saxReader) {
this.is2D = true;
if (parent == null) {
this.processXml2 (this, saxReader);
parent = this;
} else {
this.no3D = parent.checkFilterKey ("NO3D");
this.noHydrogens = parent.noHydrogens;
this.processXml2 (parent, saxReader);
this.filter = parent.filter;
}}, "J.adapter.readers.xml.XmlReader,~O");
Clazz.overrideMethod (c$, "processStartElement", 
function (localName, nodeName) {
var id = this.atts.get ("id");
if ("fragment".equals (localName)) {
this.fragments.push (this.thisFragment = id);
return;
}if ("n".equals (localName)) {
this.setNode (id);
return;
}if ("t".equals (localName)) {
this.textBuffer = "";
}if ("s".equals (localName)) {
this.setKeepChars (true);
}if ("b".equals (localName)) {
this.setBond ();
return;
}}, "~S,~S");
Clazz.overrideMethod (c$, "processEndElement", 
function (localName) {
if ("fragment".equals (localName)) {
this.thisFragment = this.fragments.pop ();
return;
}if ("n".equals (localName)) {
this.thisNode = (this.nodes.size () == 0 ? null : this.nodes.pop ());
return;
}if ("s".equals (localName)) {
this.textBuffer += this.chars.toString ();
}if ("t".equals (localName)) {
if (this.thisNode == null) {
System.out.println ("XmlChemDrawReader unassigned text: " + this.textBuffer);
} else {
this.thisNode.text = this.textBuffer;
if (this.atom.elementNumber == 0) {
System.err.println ("XmlChemDrawReader: Problem with \"" + this.textBuffer + "\"");
}if (this.thisNode.warning != null) this.parent.appendLoadNote ("Warning: " + this.textBuffer + " " + this.thisNode.warning);
}this.textBuffer = "";
}this.setKeepChars (false);
}, "~S");
Clazz.defineMethod (c$, "setNode", 
 function (id) {
var nodeType = this.atts.get ("nodetype");
if (this.asc.bsAtoms == null) this.asc.bsAtoms =  new JU.BS ();
if (this.thisNode != null) this.nodes.push (this.thisNode);
if ("_".equals (nodeType)) {
this.atom = this.thisNode = null;
return;
}this.atom = this.thisNode = Clazz.innerTypeInstance (J.adapter.readers.xml.XmlChemDrawReader.CDNode, this, null, id, nodeType, this.thisFragment, this.thisNode);
this.asc.addAtomWithMappedSerialNumber (this.atom);
this.asc.bsAtoms.set (this.atom.index);
var w = this.atts.get ("warning");
if (w != null) {
this.thisNode.warning = JU.PT.rep (w, "&apos;", "'");
this.thisNode.isValid = (w.indexOf ("ChemDraw can't interpret") < 0);
}var element = this.atts.get ("element");
var s = this.atts.get ("genericnickname");
if (s != null) {
element = s;
}this.atom.elementNumber = (!this.checkWarningOK (w) ? 0 : element == null ? 6 : this.parseIntStr (element));
element = J.api.JmolAdapter.getElementSymbol (this.atom.elementNumber);
s = this.atts.get ("isotope");
if (s != null) element = s + element;
this.setElementAndIsotope (this.atom, element);
s = this.atts.get ("charge");
if (s != null) {
this.atom.formalCharge = this.parseIntStr (s);
}var hasXYZ = (this.atts.containsKey ("xyz"));
var hasXY = (this.atts.containsKey ("p"));
if (hasXYZ && (!this.no3D || !hasXY)) {
this.is2D = false;
this.setAtom ("xyz");
} else if (this.atts.containsKey ("p")) {
this.setAtom ("p");
}s = this.atts.get ("attachments");
if (s != null) {
this.thisNode.setMultipleAttachments (JU.PT.split (s.trim (), " "));
}if (JU.Logger.debugging) JU.Logger.info ("XmlChemDraw id=" + id + " " + element + " " + this.atom);
}, "~S");
Clazz.defineMethod (c$, "checkWarningOK", 
 function (warning) {
return (warning == null || warning.indexOf ("valence") >= 0 || warning.indexOf ("very close") >= 0 || warning.indexOf ("two identical colinear bonds") >= 0);
}, "~S");
Clazz.defineMethod (c$, "setBond", 
 function () {
var atom1 = this.atts.get ("b");
var atom2 = this.atts.get ("e");
var a = this.atts.get ("beginattach");
var beginAttach = (a == null ? 0 : this.parseIntStr (a));
a = this.atts.get ("endattach");
var endAttach = (a == null ? 0 : this.parseIntStr (a));
var s = this.atts.get ("order");
var disp = this.atts.get ("display");
var disp2 = this.atts.get ("display2");
var order = 131071;
var invertEnds = false;
if (disp == null) {
if (s == null) {
order = 1;
} else {
if (s.indexOf (".") > 0 && !"Dash".equals (disp2)) {
s = s.substring (0, s.indexOf ("."));
}order = JU.Edge.getBondOrderFromString (s);
}} else if (disp.equals ("WedgeBegin")) {
order = 1025;
} else if (disp.equals ("Hash") || disp.equals ("WedgedHashBegin")) {
order = 1041;
} else if (disp.equals ("WedgeEnd")) {
invertEnds = true;
order = 1025;
} else if (disp.equals ("WedgedHashEnd")) {
invertEnds = true;
order = 1041;
} else if (disp.equals ("Wavy")) {
order = 1057;
}if (order == 131071) {
System.err.println ("XmlChemDrawReader ignoring bond type " + s);
return;
}var b = (invertEnds ? Clazz.innerTypeInstance (J.adapter.readers.xml.XmlChemDrawReader.CDBond, this, null, atom2, atom1, order) : Clazz.innerTypeInstance (J.adapter.readers.xml.XmlChemDrawReader.CDBond, this, null, atom1, atom2, order));
var node1 = this.asc.atoms[b.atomIndex1];
var node2 = this.asc.atoms[b.atomIndex2];
if (order == 1057) {
if (!this.nostereo.contains (node1)) this.nostereo.add (node1);
if (!this.nostereo.contains (node2)) this.nostereo.add (node2);
}if (node1.hasMultipleAttachments) {
node1.attachedAtom = node2;
return;
} else if (node2.hasMultipleAttachments) {
node2.attachedAtom = node1;
return;
}if (node1.isFragment && beginAttach == 0) beginAttach = 1;
if (node2.isFragment && endAttach == 0) endAttach = 1;
if (beginAttach > 0) {
(invertEnds ? node2 : node1).addAttachedAtom (b, beginAttach);
}if (endAttach > 0) {
(invertEnds ? node1 : node2).addAttachedAtom (b, endAttach);
}if (node1.isExternalPt) {
node1.internalBond = b;
node2.parentNode.addExternalPoint (node2);
}if (node2.isExternalPt) {
node2.internalBond = b;
node1.parentNode.addExternalPoint (node1);
}this.asc.addBondNoCheck (b);
});
Clazz.defineMethod (c$, "setAtom", 
 function (key) {
var xyz = this.atts.get (key);
var tokens = JU.PT.getTokens (xyz);
var x = this.parseFloatStr (tokens[0]);
var y = -this.parseFloatStr (tokens[1]);
var z = (key === "xyz" ? this.parseFloatStr (tokens[2]) : 0);
if (x < this.minX) this.minX = x;
if (x > this.maxX) this.maxX = x;
if (y < this.minY) this.minY = y;
if (y > this.maxY) this.maxY = y;
if (z < this.minZ) this.minZ = z;
if (z > this.maxZ) this.maxZ = z;
this.atom.set (x, y, z);
}, "~S");
Clazz.overrideMethod (c$, "finalizeSubclassReader", 
function () {
this.fixConnections ();
this.fixInvalidAtoms ();
this.centerAndScale ();
this.asc.setInfo ("minimize3D", Boolean.$valueOf (!this.is2D && !this.noHydrogens));
this.set2D ();
this.asc.setInfo ("is2D", Boolean.$valueOf (this.is2D));
if (this.is2D) this.asc.setModelInfoForSet ("dimension", "2D", this.asc.iSet);
this.parent.appendLoadNote ("ChemDraw CDXML: " + (this.is2D ? "2D" : "3D"));
});
Clazz.defineMethod (c$, "fixConnections", 
 function () {
for (var i = this.asc.ac; --i >= 0; ) {
var a = this.asc.atoms[i];
if (a.isFragment || a.hasMultipleAttachments) a.fixAttachments ();
}
for (var i = 0, n = this.asc.bondCount; i < n; i++) {
var b = this.asc.bonds[i];
if (b == null) {
continue;
}var a1 = this.asc.atoms[b.atomIndex1];
var a2 = this.asc.atoms[b.atomIndex2];
a1.isConnected = true;
a2.isConnected = true;
if (this.nostereo.contains (a1) != this.nostereo.contains (a2)) {
b.order = 1;
}}
});
Clazz.defineMethod (c$, "centerAndScale", 
 function () {
if (this.minX > this.maxX) return;
var sum = 0;
var n = 0;
var lenH = 1;
for (var i = this.asc.bondCount; --i >= 0; ) {
var a1 = this.asc.atoms[this.asc.bonds[i].atomIndex1];
var a2 = this.asc.atoms[this.asc.bonds[i].atomIndex2];
var d = a1.distance (a2);
if (a1.elementNumber > 1 && a2.elementNumber > 1) {
sum += d;
n++;
} else {
lenH = d;
}}
var f = (sum > 0 ? 1.45 * n / sum : lenH > 0 ? 1 / lenH : 1);
if (f > 0.5) f = 1;
var cx = (this.maxX + this.minX) / 2;
var cy = (this.maxY + this.minY) / 2;
var cz = (this.maxZ + this.minZ) / 2;
for (var i = this.asc.ac; --i >= 0; ) {
var a = this.asc.atoms[i];
a.x = (a.x - cx) * f;
a.y = (a.y - cy) * f;
a.z = (a.z - cz) * f;
}
});
Clazz.defineMethod (c$, "fixInvalidAtoms", 
 function () {
for (var i = this.asc.ac; --i >= 0; ) {
var a = this.asc.atoms[i];
a.atomSerial = -2147483648;
if (a.isFragment || a.isExternalPt || !a.isConnected && (!a.isValid || a.elementNumber == 6 || a.elementNumber == 0)) {
this.asc.bsAtoms.clear (a.index);
}}
});
c$.$XmlChemDrawReader$CDNode$ = function () {
Clazz.pu$h(self.c$);
c$ = Clazz.decorateAsClass (function () {
Clazz.prepareCallback (this, arguments);
this.warning = null;
this.id = null;
this.intID = 0;
this.isValid = true;
this.isConnected = false;
this.isFragment = false;
this.isExternalPt = false;
this.nodeType = null;
this.fragment = null;
this.text = null;
this.parentNode = null;
this.orderedExternalPoints = null;
this.orderedAttachedBonds = null;
this.internalBond = null;
this.attachments = null;
this.hasMultipleAttachments = false;
this.attachedAtom = null;
this.isGeneric = false;
Clazz.instantialize (this, arguments);
}, J.adapter.readers.xml.XmlChemDrawReader, "CDNode", J.adapter.smarter.Atom);
Clazz.makeConstructor (c$, 
function (a, b, c, d) {
Clazz.superConstructor (this, J.adapter.readers.xml.XmlChemDrawReader.CDNode, []);
this.id = a;
this.fragment = c;
this.atomSerial = this.intID = Integer.parseInt (a);
this.nodeType = b;
this.parentNode = d;
this.isFragment = "Fragment".equals (b) || "Nickname".equals (b);
this.isExternalPt = "ExternalConnectionPoint".equals (b);
this.isGeneric = "GenericNickname".equals (b);
}, "~S,~S,~S,J.adapter.readers.xml.XmlChemDrawReader.CDNode");
Clazz.defineMethod (c$, "setMultipleAttachments", 
function (a) {
this.attachments = a;
this.hasMultipleAttachments = true;
}, "~A");
Clazz.defineMethod (c$, "addExternalPoint", 
function (a) {
if (this.orderedExternalPoints == null) this.orderedExternalPoints =  new JU.Lst ();
var b = this.orderedExternalPoints.size ();
while (--b >= 0 && this.orderedExternalPoints.get (b).intID >= a.intID) {
}
this.orderedExternalPoints.add (++b, a);
}, "J.adapter.readers.xml.XmlChemDrawReader.CDNode");
Clazz.defineMethod (c$, "addAttachedAtom", 
function (a, b) {
if (this.orderedAttachedBonds == null) this.orderedAttachedBonds =  new JU.Lst ();
var c = this.orderedAttachedBonds.size ();
while (--c >= 0 && (this.orderedAttachedBonds.get (c)[0]).intValue () > b) {
}
this.orderedAttachedBonds.add (++c,  Clazz.newArray (-1, [Integer.$valueOf (b), a]));
}, "J.adapter.readers.xml.XmlChemDrawReader.CDBond,~N");
Clazz.defineMethod (c$, "fixAttachments", 
function () {
if (this.hasMultipleAttachments && this.attachedAtom != null) {
var a = JU.Edge.getBondOrderFromString ("partial");
var b = this.attachedAtom.index;
for (var c = this.attachments.length; --c >= 0; ) {
var d = this.b$["J.adapter.readers.xml.XmlChemDrawReader"].asc.getAtomFromName (this.attachments[c]);
if (d != null) this.b$["J.adapter.readers.xml.XmlChemDrawReader"].asc.addBondNoCheck ( new J.adapter.smarter.Bond (b, d.index, a));
}
}if (this.orderedExternalPoints == null || this.text == null) return;
var a = this.orderedExternalPoints.size ();
if (a != this.orderedAttachedBonds.size ()) {
System.err.println ("cannot fix attachments for " + this.text);
}for (var b = a; --b >= 0; ) {
var c = this.orderedExternalPoints.get (b);
var d = this.orderedAttachedBonds.get (b)[1];
var e = (d.atomIndex2 == this.index);
var f = d.getOtherNode (this);
var g = c;
this.b$["J.adapter.readers.xml.XmlChemDrawReader"].asc.addBondNoCheck ( new J.adapter.smarter.Bond ((e ? f.index : g.index), (e ? g.index : f.index), d.order));
}
});
Clazz.overrideMethod (c$, "toString", 
function () {
return this.id + " " + this.elementSymbol + " " + this.elementNumber + " index=" + this.index + " ext=" + this.isExternalPt + " frag=" + this.isFragment + " " + this.elementSymbol + " " + this.x + " " + this.y;
});
c$ = Clazz.p0p ();
};
c$.$XmlChemDrawReader$CDBond$ = function () {
Clazz.pu$h(self.c$);
c$ = Clazz.decorateAsClass (function () {
Clazz.prepareCallback (this, arguments);
this.id1 = null;
this.id2 = null;
Clazz.instantialize (this, arguments);
}, J.adapter.readers.xml.XmlChemDrawReader, "CDBond", J.adapter.smarter.Bond);
Clazz.makeConstructor (c$, 
function (a, b, c) {
Clazz.superConstructor (this, J.adapter.readers.xml.XmlChemDrawReader.CDBond, [this.b$["J.adapter.readers.xml.XmlChemDrawReader"].asc.getAtomFromName (a).index, this.b$["J.adapter.readers.xml.XmlChemDrawReader"].asc.getAtomFromName (b).index, c]);
this.id1 = a;
this.id2 = b;
}, "~S,~S,~N");
Clazz.defineMethod (c$, "getOtherNode", 
function (a) {
return this.b$["J.adapter.readers.xml.XmlChemDrawReader"].asc.atoms[this.atomIndex1 == a.index ? this.atomIndex2 : this.atomIndex1];
}, "J.adapter.readers.xml.XmlChemDrawReader.CDNode");
Clazz.defineMethod (c$, "toString", 
function () {
return Clazz.superCall (this, J.adapter.readers.xml.XmlChemDrawReader.CDBond, "toString", []) + " " + this.id1 + " " + this.id2;
});
c$ = Clazz.p0p ();
};
});
