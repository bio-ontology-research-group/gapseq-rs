//! XML namespace URIs and prefixes used in the emitted document.

pub const CORE: &str = "http://www.sbml.org/sbml/level3/version1/core";
pub const FBC: &str = "http://www.sbml.org/sbml/level3/version1/fbc/version2";
pub const GROUPS: &str = "http://www.sbml.org/sbml/level3/version1/groups/version1";
pub const XHTML: &str = "http://www.w3.org/1999/xhtml";

#[allow(dead_code)]
pub const FBC_PREFIX: &str = "fbc";
#[allow(dead_code)]
pub const GROUPS_PREFIX: &str = "groups";

// Pre-built prefixed tag / attribute names. quick-xml's `BytesStart` borrows
// its name; using `&'static str` avoids per-event `format!` allocations
// (and sidesteps lifetime headaches).
pub const ATTR_XMLNS_FBC: &str = "xmlns:fbc";
pub const ATTR_XMLNS_GROUPS: &str = "xmlns:groups";
pub const ATTR_FBC_REQUIRED: &str = "fbc:required";
pub const ATTR_GROUPS_REQUIRED: &str = "groups:required";
pub const ATTR_FBC_STRICT: &str = "fbc:strict";
pub const ATTR_FBC_CHARGE: &str = "fbc:charge";
pub const ATTR_FBC_FORMULA: &str = "fbc:chemicalFormula";
pub const ATTR_FBC_LOWER: &str = "fbc:lowerFluxBound";
pub const ATTR_FBC_UPPER: &str = "fbc:upperFluxBound";
pub const ATTR_FBC_GENEPRODUCT: &str = "fbc:geneProduct";
pub const ATTR_FBC_ID: &str = "fbc:id";
pub const ATTR_FBC_TYPE: &str = "fbc:type";
pub const ATTR_FBC_REACTION: &str = "fbc:reaction";
pub const ATTR_FBC_COEFFICIENT: &str = "fbc:coefficient";
pub const ATTR_FBC_LABEL: &str = "fbc:label";
pub const ATTR_FBC_ACTIVE_OBJECTIVE: &str = "fbc:activeObjective";
pub const ATTR_GROUPS_ID: &str = "groups:id";
pub const ATTR_GROUPS_NAME: &str = "groups:name";
pub const ATTR_GROUPS_KIND: &str = "groups:kind";
pub const ATTR_GROUPS_IDREF: &str = "groups:idRef";

pub const TAG_FBC_GPA: &str = "fbc:geneProductAssociation";
pub const TAG_FBC_GENEPRODUCTREF: &str = "fbc:geneProductRef";
pub const TAG_FBC_AND: &str = "fbc:and";
pub const TAG_FBC_OR: &str = "fbc:or";
pub const TAG_FBC_OBJECTIVES: &str = "fbc:listOfObjectives";
pub const TAG_FBC_OBJECTIVE: &str = "fbc:objective";
pub const TAG_FBC_FLUX_OBJECTIVES: &str = "fbc:listOfFluxObjectives";
pub const TAG_FBC_FLUX_OBJECTIVE: &str = "fbc:fluxObjective";
pub const TAG_FBC_GENEPRODUCTS: &str = "fbc:listOfGeneProducts";
pub const TAG_FBC_GENEPRODUCT: &str = "fbc:geneProduct";
pub const TAG_GROUPS_GROUPS: &str = "groups:listOfGroups";
pub const TAG_GROUPS_GROUP: &str = "groups:group";
pub const TAG_GROUPS_MEMBERS: &str = "groups:listOfMembers";
pub const TAG_GROUPS_MEMBER: &str = "groups:member";

/// SBO term for a flux balance constraint parameter.
pub const SBO_FLUX_BOUND: &str = "SBO:0000626";
/// SBO term for default flux bound.
#[allow(dead_code)]
pub const SBO_DEFAULT_FLUX_BOUND: &str = "SBO:0000626";
