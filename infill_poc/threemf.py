"""Writer for PrusaSlicer-flavoured 3MF files with per-volume setting overrides.

Format reverse-engineered from a PrusaSlicer 2.5.0 `--export-3mf`: every volume of an object is
concatenated into a single mesh in 3D/3dmodel.model, and Metadata/Slic3r_PE_model.config assigns
triangle index ranges (firstid..lastid) to volumes. On import, any volume metadata key that is not
one of PrusaSlicer's reserved keys is applied as a print setting for that volume, which is how a
`fill_density` override reaches the slicer. Bambu Studio / OrcaSlicer use a different file
(Metadata/model_settings.config) and are not targeted here.
"""
import xml.etree.ElementTree as ET
import zipfile
from xml.sax.saxutils import escape, quoteattr

import numpy as np

IDENTITY_4x4 = "1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1"

CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
 <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
 <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>
</Types>
"""

RELS = """<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
 <Relationship Target="/3D/3dmodel.model" Id="rel-1" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>
</Relationships>
"""


def _fmt(v):
    return f"{v:.6f}".rstrip("0").rstrip(".") or "0"


def write_prusa_3mf(path, base_mesh, modifiers, translation, object_name="probe_part"):
    """base_mesh: trimesh part. modifiers: list of (name, trimesh, {setting_key: value}).
    translation: (x, y, z) placement of the object on the bed."""
    volumes = [("base", base_mesh, "ModelPart", {})]
    volumes += [(name, mesh, "ParameterModifier", settings) for name, mesh, settings in modifiers]

    vert_xml, tri_xml, vol_xml = [], [], []
    v_offset = t_offset = 0
    for name, mesh, vtype, settings in volumes:
        verts = np.asarray(mesh.vertices, dtype=float)
        faces = np.asarray(mesh.faces, dtype=int) + v_offset
        vert_xml.extend(f'     <vertex x="{_fmt(x)}" y="{_fmt(y)}" z="{_fmt(z)}"/>' for x, y, z in verts)
        tri_xml.extend(f'     <triangle v1="{a}" v2="{b}" v3="{c}"/>' for a, b, c in faces)

        meta = [("name", name), ("volume_type", vtype), ("matrix", IDENTITY_4x4)]
        if vtype == "ParameterModifier":
            meta.insert(1, ("modifier", "1"))
        meta += list(settings.items())
        lines = [f'  <volume firstid="{t_offset}" lastid="{t_offset + len(faces) - 1}">']
        lines += [f'   <metadata type="volume" key={quoteattr(k)} value={quoteattr(str(v))}/>' for k, v in meta]
        lines += ['   <mesh edges_fixed="0" degenerate_facets="0" facets_removed="0" facets_reversed="0" backwards_edges="0"/>',
                  "  </volume>"]
        vol_xml.extend(lines)
        v_offset += len(verts)
        t_offset += len(faces)

    tx, ty, tz = translation
    model = "\n".join([
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<model unit="millimeter" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" '
        'xmlns:slic3rpe="http://schemas.slic3r.org/3mf/2017/06">',
        ' <metadata name="slic3rpe:Version3mf">1</metadata>',
        f' <metadata name="Title">{escape(object_name)}</metadata>',
        ' <metadata name="Application">stiffness_infill_poc</metadata>',
        ' <resources>',
        '  <object id="1" type="model">',
        '   <mesh>',
        '    <vertices>', *vert_xml, '    </vertices>',
        '    <triangles>', *tri_xml, '    </triangles>',
        '   </mesh>',
        '  </object>',
        ' </resources>',
        ' <build>',
        f'  <item objectid="1" transform="1 0 0 0 1 0 0 0 1 {_fmt(tx)} {_fmt(ty)} {_fmt(tz)}" printable="1"/>',
        ' </build>',
        '</model>',
        '',
    ])
    config = "\n".join([
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<config>',
        ' <object id="1" instances_count="1">',
        f'  <metadata type="object" key="name" value={quoteattr(object_name)}/>',
        *vol_xml,
        ' </object>',
        '</config>',
        '',
    ])
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", CONTENT_TYPES)
        z.writestr("_rels/.rels", RELS)
        z.writestr("3D/3dmodel.model", model)
        z.writestr("Metadata/Slic3r_PE_model.config", config)


def read_volume_settings(path):
    """Parse Slic3r_PE_model.config -> list of {key: value} per volume (used for round-trip checks)."""
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("Metadata/Slic3r_PE_model.config"))
    return [{m.get("key"): m.get("value") for m in vol.findall("metadata")} for vol in root.iter("volume")]
