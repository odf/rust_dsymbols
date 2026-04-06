use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use cgmath::prelude::*;
use cgmath::{point3, vec3, Point3};
use three_d::{Mat4, Vec3};

use rust_dsymbols::delaney3d::pseudo_toroidal_cover;
use rust_dsymbols::display::mesh::{ItemType, Mesh, decompose_mesh, scaled_mesh};
use rust_dsymbols::dsets::DSet;
use rust_dsymbols::dsyms::{DSym, PartialDSym};
use rust_dsymbols::geometry::vec_matrix::VecMatrix;
use rust_dsymbols::tilings::{
    Skeleton, chamber_positions, gram_matrix, invariant_basis, tile_surfaces
};


#[derive(Clone, Copy, PartialEq)]
struct Options {
    tile_scale: f64,
    edge_radius: f64,
    vertex_color: three_d::egui::Color32,
    edge_color: three_d::egui::Color32,
    face_color: three_d::egui::Color32,
    background_color: three_d::egui::Color32,
    sun_direction: Vec3,
    sun_casts_shadows: bool,
}


impl Default for Options {
    fn default() -> Self {
        Self {
            tile_scale: 0.75,
            edge_radius: 0.05,
            vertex_color: three_d::egui::Color32::BLACK,
            edge_color: three_d::egui::Color32::BLUE,
            face_color: three_d::egui::Color32::RED,
            background_color: three_d::egui::Color32::GRAY,
            sun_direction: vec3(1.0, -1.0, -1.0),
            sun_casts_shadows: false,
        }
    }
}


pub struct Tiling {
    spec: String,
    ds: RefCell<Option<Rc<Result<PartialDSym, String>>>>,
    cov: RefCell<Option<Rc<Result<PartialDSym, String>>>>,
    skel: RefCell<Option<Rc<Result<Skeleton, String>>>>,
}


impl Tiling {
    pub fn from(spec: &str) -> Tiling {
        Tiling {
            spec: String::from(spec),
            ds: RefCell::new(None),
            cov: RefCell::new(None),
            skel: RefCell::new(None)
        }
    }

    pub fn ds(&self) -> Rc<Result<PartialDSym, String>> {
        if self.ds.borrow().is_none() {
            let maybe_ds = self.spec.parse::<PartialDSym>();
            self.ds.replace(Some(Rc::new(maybe_ds)));
        }

        self.ds.borrow().as_ref().unwrap().clone()
    }

    pub fn cov(&self) -> Rc<Result<PartialDSym, String>> {
        if self.cov.borrow().is_none() {
            let maybe_cov = self.cov_impl();
            self.cov.replace(Some(Rc::new(maybe_cov)));
        }

        self.cov.borrow().as_ref().unwrap().clone()
    }

    pub fn skel(&self) -> Rc<Result<Skeleton, String>> {
        if self.skel.borrow().is_none() {
            let maybe_skel = match self.cov().as_ref() {
                Ok(val) => Ok(Skeleton::of(val)),
                Err(msg) => Err(msg.clone()),
            };
            self.skel.replace(Some(Rc::new(maybe_skel)));
        }

        self.skel.borrow().as_ref().unwrap().clone()
    }

    pub fn clear_cache(&self) {
        self.ds.replace(None);
        self.cov.replace(None);
        self.skel.replace(None);
    }

    fn cov_impl(&self) -> Result<PartialDSym, String> {
        match self.ds().as_ref() {
            Ok(ds) => {
                if ds.dim() != 3 {
                    Err("is not three_dimensional".to_string())
                } else if !ds.is_complete() {
                    Err("is incomplete".to_string())
                } else {
                    for i in 0..ds.dim() {
                        for d in ds.orbit_reps_2d(i, i + 1) {
                            let v = ds.v(i, i + 1, d).unwrap();
                            if v > 6 && v == 5 {
                                return Err(String::from(
                                    "violates the crystallographic restriction"
                                ));
                            }
                        }
                    }

                    pseudo_toroidal_cover(ds).ok_or(
                        "does not have a pseudo-toroidal cover".to_string()
                    )
                }
            },
            Err(msg) => Err(msg.clone())
        }
    }
}


struct State {
    options: Options,
    previous_options: Options,
    catalog: HashMap<String, Vec<Tiling>>,
    collection_name: String,
    index_in_collection: usize,
    camera: three_d::Camera,
    models: Vec<three_d::Gm<three_d::InstancedMesh, three_d::PhysicalMaterial>>,
}


fn main() {
    // On the web, this creates a canvas instead.
    let window = three_d::Window::new(three_d::WindowSettings {
        title: "Rust 3d Test".to_string(),
        max_size: Some((1280, 720)),
        ..Default::default()
    })
    .unwrap();

    let mut context = window.gl();
    let mut gui = three_d::GUI::new(&context);

    let camera = three_d::Camera::new_perspective(
        window.viewport(),
        vec3(0.0, 2.0, 8.0),
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 1.0, 0.0),
        three_d::degrees(25.0),
        0.1,
        1000.0,
    );

    let options = Default::default();
    let previous_options = options;

    let builtin = [
        /* bcu */ "<1.1:2 3:2,1 2,1 2,2:4,4 2,6>",
        /* pcu */ "<1.1:1 3:1,1,1,1:4,3,4>",
        /* nbo */ "<1.1:2 3:2,1 2,1 2,2:6,4 2,4>",
        /* dia */ "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>",
        /* srs */ "<1.1:10 3:2 4 6 8 10,10 3 5 7 9,10 9 8 7 6,4 3 10 9 8:10,3 2 2,10>"
    ];

    let catalog = HashMap::from([
        (
            String::from("__builtin__"),
            builtin.into_iter().map(Tiling::from).collect::<Vec<_>>()
        )
    ]);

    let collection_name = "__builtin__".to_string();
    let index_in_collection = 4;

    #[cfg(feature = "pprof")]
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build().unwrap();

    let models = build_models(
        &context, &catalog[&collection_name[..]][index_in_collection], &options
    );

    #[cfg(feature = "pprof")]
    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create("flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
    };

    let mut state = State {
        options, previous_options,
        catalog, collection_name, index_in_collection,
        camera, models,
    };

    window.render_loop(move |mut frame_input| {
        render_callback(&mut frame_input, &mut context, &mut gui, &mut state)
    });
}


fn render_callback(
    frame_input: &mut three_d::FrameInput,
    context: &mut three_d::Context,
    gui: &mut three_d::GUI,
    state: &mut State,
)
    -> three_d::FrameOutput
{
    let mut panel_width = 0.0;

    gui.update(
        &mut frame_input.events,
        frame_input.accumulated_time,
        frame_input.viewport,
        frame_input.device_pixel_ratio,
        |gui_context| {
            panel_width = gui_callback(state, gui_context);
        },
    );

    let wd = panel_width * frame_input.device_pixel_ratio;

    let viewport = three_d::Viewport {
        x: wd as i32,
        y: 0,
        width: frame_input.viewport.width - wd as u32,
        height: frame_input.viewport.height,
    };
    state.camera.set_viewport(viewport);

    rust_dsymbols::display::controls::orbit_control_update_camera(
        &mut state.camera, &mut frame_input.events, 1.0, 1000.0
    );

    let white = three_d::Srgba::WHITE;
    let sun_dir = (
        state.camera.view().invert().unwrap() *
        state.options.sun_direction.extend(0.0)
    ).truncate();

    let mut sun = three_d::DirectionalLight::new(context, 2.0, white, sun_dir);
    let ambient = three_d::AmbientLight::new(context, 0.1, white);

    if state.options != state.previous_options {
        state.models = build_models(
            context,
            &state.catalog[&state.collection_name][state.index_in_collection],
            &state.options
        );
        state.previous_options = state.options;
    }

    if state.options.sun_casts_shadows {
        sun.generate_shadow_map(9192, &state.models);
    }

    let [r, g, b, a] = state.options.background_color.to_normalized_gamma_f32();

    frame_input.screen()
        .clear(three_d::ClearState::color_and_depth(r, g, b, a, 1.0))
        .render(&state.camera, &state.models, &[&sun, &ambient])
        .write(|| gui.render()).unwrap();

    // Ensures a valid return value.
    three_d::FrameOutput::default()
}


fn gui_callback(state: &mut State, gui_context: &three_d::egui::Context)
    -> f32
{
    three_d::egui::SidePanel::left("side_panel").show(gui_context, |ui| {
        let options = &mut state.options;

        ui.heading("Tiling");
        ui.add_space(12.0);

        ui.label(format!("Collection '{}'", state.collection_name));
        ui.label(format!("Index {}", state.index_in_collection));
        ui.add_space(24.0);

        ui.heading("Appearance");
        ui.add_space(12.0);

        ui.add(
            three_d::egui::Slider::new(&mut options.tile_scale, 0.1..=1.0)
                .text("Tile scale")
        );
        ui.add(
            three_d::egui::Slider::new(&mut options.edge_radius, 0.0..=0.1)
                .text("Edge radius")
        );
        ui.add_space(12.0);

        ui.horizontal(|ui| {
            ui.color_edit_button_srgba(&mut options.face_color);
            ui.label("Face color");
        });
        ui.horizontal(|ui| {
            ui.color_edit_button_srgba(&mut options.edge_color);
            ui.label("Edge color");
        });
        ui.horizontal(|ui| {
            ui.color_edit_button_srgba(&mut options.vertex_color);
            ui.label("Vertex color");
        });
        ui.add_space(12.0);

        ui.checkbox(&mut options.sun_casts_shadows, "Shadows On");
        ui.add_space(12.0);

        ui.horizontal(|ui| {
            ui.color_edit_button_srgba(&mut options.background_color);
            ui.label("Background color");
        });
    });
    gui_context.used_rect().width()
}


fn build_models(context: &three_d::Context, til: &Tiling, options: &Options)
    -> Vec<three_d::Gm<three_d::InstancedMesh, three_d::PhysicalMaterial>>
{
    let instances = three_d::Instances {
        transformations: vec![
            Mat4::from_scale(1.0),
        ],
        ..Default::default()
    };

    tiles(til).iter().flat_map(|tile_mesh|
        decompose_mesh(
                &scaled_mesh(tile_mesh, options.tile_scale),
                options.edge_radius
            )
            .iter()
            .map(|(part_mesh, item_type)| {
                let color = match item_type {
                    ItemType::Vertex => options.vertex_color,
                    ItemType::Edge => options.edge_color,
                    ItemType::Face => options.face_color,
                };

                three_d::Gm::new(
                    three_d::InstancedMesh::new(
                        context,
                        &instances,
                        &part_mesh.to_cpu_mesh()
                    ),
                    three_d::PhysicalMaterial {
                        albedo: color.to_normalized_gamma_f32().into(),
                        metallic: 0.0,
                        roughness: 0.5,
                        ..Default::default()
                    }
                )
            })
            .collect::<Vec<_>>()
        )
        .collect()
}


fn tiles(til: &Tiling) -> Vec<Mesh<Point3<f64>>> {
    let ds = til.ds();
    let ds = ds.as_ref().as_ref().unwrap();
    let cov = til.cov();
    let cov = cov.as_ref().as_ref().unwrap();
    let skel = til.skel();
    let skel = skel.as_ref().as_ref().unwrap();

    let basis = invariant_basis(&gram_matrix(ds, cov, skel).unwrap()
    ).transpose();
    let pos = skel.graph.vertices().iter()
        .map(|&v| (v, skel.graph.position(v)))
        .collect();

    let (scale, shift) = scale_and_shift(cov, skel, &basis);
    let reps = cov.orbit_reps([0, 1, 2], cov.elements());

    tile_surfaces(cov, skel, &pos, reps).iter().map(|(vertices, faces)| {
        let vs = vertices.iter().map(|v| {
            let v = scale * &basis * v.to_f64().unwrap() + &shift;
            point3(v[(0, 0)], v[(1, 0)], v[(2, 0)])
        });

        Mesh::from_oriented_faces(vs, faces.clone()).unwrap()
    }).collect()
}


fn scale_and_shift(cov: &PartialDSym, skel: &Skeleton, basis: &VecMatrix<f64>)
    -> (f64, VecMatrix<f64>)
{
    let ch_pos = chamber_positions(cov, skel);
    let mut s = 0.0;
    let mut n = 0;

    for d in  cov.orbit_reps([0, 2, 3], cov.elements()) {
        let v = 2.0 * basis * (&ch_pos[&(d, 1)] - &ch_pos[&(d, 0)]).to_f64().unwrap();
        s += v.norm();
        n += 1;
    };

    let scale = n as f64 / s;
    let shift = -scale * basis * ch_pos[&(1, 3)].to_f64().unwrap();

    (scale, shift)
}


trait ToCpuMesh {
    fn to_cpu_mesh(&self) -> three_d::CpuMesh;
}


impl ToCpuMesh for Mesh<Point3<f64>> {
    fn to_cpu_mesh(&self) -> three_d::CpuMesh {
        let trimesh = self.triangulate().unwrap();

        let positions: Vec<_> = trimesh.vertices().iter()
            .map(|&p| p.to_vec())
            .collect();

        let indices: Vec<_> = trimesh.face_indices().iter()
            .flatten()
            .map(|&x| x as u32)
            .collect();

        let mut mesh = three_d::CpuMesh {
            positions: three_d::Positions::F64(positions),
            indices: three_d::Indices::U32(indices),
            ..Default::default()
        };
        mesh.compute_normals();
        mesh
    }
}
