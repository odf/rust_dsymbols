use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::rc::Rc;

use cgmath::prelude::*;
use cgmath::{point3, vec3, Point3};
use egui_file::FileDialog;
use lru::LruCache;
use sha2::{Sha256, Digest};
use three_d::{Mat4, Vec3};

use rust_dsymbols::delaney3d::{
    obeys_crystallographic_restriction, pseudo_toroidal_cover
};
use rust_dsymbols::display::mesh::{ItemType, Mesh, decompose_mesh, scaled_mesh};
use rust_dsymbols::dsets::DSet;
use rust_dsymbols::dsyms::PartialDSym;
use rust_dsymbols::geometry::vec_matrix::VecMatrix;
use rust_dsymbols::tilings::{
    Skeleton, chamber_positions, gram_matrix, invariant_basis, tile_surfaces
};


pub struct Tiling {
    ds: Result<PartialDSym, String>,
    cov: RefCell<Option<Rc<Result<PartialDSym, String>>>>,
    skel: RefCell<Option<Rc<Result<Skeleton, String>>>>,
}


impl Tiling {
    pub fn from(spec: &str) -> Tiling {
        Tiling {
            ds: spec.parse::<PartialDSym>(),
            cov: RefCell::new(None),
            skel: RefCell::new(None)
        }
    }

    pub fn ds(&self) -> &Result<PartialDSym, String> {
        &self.ds
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

    fn cov_impl(&self) -> Result<PartialDSym, String> {
        let ds = self.ds().as_ref()?;

        if ds.dim() != 3 {
            Err("is not three_dimensional".to_string())
        } else if !ds.is_complete() {
            Err("is incomplete".to_string())
        } else if !obeys_crystallographic_restriction(ds) {
            Err("violates the crystallographic restriction".to_string())
        } else {
            pseudo_toroidal_cover(ds).ok_or(
                "does not have a pseudo-toroidal cover".to_string()
            )
        }
    }
}


struct Catalog<T> {
    collections: HashMap<String, Vec<T>>,
    collection_name: String,
    index_in_collection: usize,
}


impl<T> Catalog<T> {
    fn from(
        title: &str,
        collection: impl IntoIterator<Item=T>
    ) -> Self {
        let mut result = Self {
            collections: HashMap::new(),
            collection_name: Default::default(),
            index_in_collection: Default::default()
        };
        result.add_collection(title, collection);

        result
    }

    fn get(&self) -> Option<&T> {
        self.collections
            .get(&self.collection_name)?
            .get(self.index_in_collection)
    }

    fn add_collection(
        &mut self,
        title: &str,
        collection: impl IntoIterator<Item=T>
    ) {
        self.collections.insert(title.to_string(), collection.into_iter().collect());
        self.collection_name = title.to_string();
        self.index_in_collection = 0;
    }
}


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


#[derive(Eq, Hash, PartialEq)]
struct CacheKey {
    tiling_spec: String,
    tile_scale: i64,
    edge_radius: i64,
}


struct State {
    options: Options,
    catalog: Catalog<String>,
    camera: three_d::Camera,
    opened_file: Option<PathBuf>,
    open_file_dialog: Option<FileDialog>,
    tiling_cache: LruCache<String, Tiling>,
    parts_cache: LruCache<CacheKey, Result<Vec<(three_d::CpuMesh, ItemType)>, String>>,
    message: String,
}


impl State {
    fn update_caches(&mut self) {
        if let Some(tkey) = self.tiling_cache_key() {
            let spec = self.catalog.get().unwrap();
            let pkey = self.parts_cache_key().unwrap();

            if self.tiling_cache.get(&tkey).is_none() {
                let til = Tiling::from(spec);
                self.tiling_cache.put(tkey.clone(), til);
            }

            if self.parts_cache.get(&pkey).is_none() {
                let parts = build_parts(
                    self.tiling_cache.get(&tkey).unwrap(),
                    &self.options
                );
                self.parts_cache.put(pkey, parts);
            }
        }
    }

    fn tiling_cache_key(&self) -> Option<String> {
        self.catalog.get().map(|spec| hex::encode(Sha256::digest(spec)))
    }

    fn parts_cache_key(&self) -> Option<CacheKey> {
        self.catalog.get().map(|spec|
            CacheKey {
                tiling_spec: hex::encode(Sha256::digest(spec)),
                tile_scale: (self.options.tile_scale * 100.0) as i64,
                edge_radius: (self.options.edge_radius * 1000.0) as i64,
            }
        )
    }
}


fn main() {
    // On the web, this creates a canvas instead.
    let window = three_d::Window::new(three_d::WindowSettings {
        title: "Rust 3d Test".to_string(),
        max_size: Some((1280, 720)),
        ..Default::default()
    })
    .unwrap();

    let camera = three_d::Camera::new_perspective(
        window.viewport(),
        vec3(0.0, 2.0, 8.0),
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 1.0, 0.0),
        three_d::degrees(25.0),
        0.1,
        1000.0,
    );

    let builtin = [
        /* bcu */ "<1.1:2 3:2,1 2,1 2,2:4,4 2,6>",
        /* pcu */ "<1.1:1 3:1,1,1,1:4,3,4>",
        /* nbo */ "<1.1:2 3:2,1 2,1 2,2:6,4 2,4>",
        /* dia */ "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>",
        /* srs */ "<1.1:10 3:2 4 6 8 10,10 3 5 7 9,10 9 8 7 6,4 3 10 9 8:10,3 2 2,10>"
    ];

    let mut state = State {
        options: Default::default(),
        catalog: Catalog::from("__builtin__", builtin.into_iter().map(String::from)),
        camera,
        opened_file: None,
        open_file_dialog: None,
        tiling_cache: LruCache::new(NonZero::new(10).unwrap()),
        parts_cache: LruCache::new(NonZero::new(10).unwrap()),
        message: "Initializing...".to_string(),
    };

    let mut context = window.gl();
    let mut gui = three_d::GUI::new(&context);

    #[cfg(feature = "pprof")] {
        profile_build_models(&mut context, &state);
    }

    window.render_loop(move |mut frame_input| {
        render_callback(&mut frame_input, &mut context, &mut gui, &mut state)
    });
}


#[cfg(feature = "pprof")]
fn profile_build_models(context: &mut three_d::Context, state: &State) {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build().unwrap();

    let tiling = Tiling::from(state.catalog.get());
    let options = &state.options;

    if let Ok(parts) = build_parts(&tiling, options) {
        build_models(&context, &parts, options);
    }

    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create("flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
    };
}


fn render_callback(
    frame_input: &mut three_d::FrameInput,
    context: &mut three_d::Context,
    gui: &mut three_d::GUI,
    state: &mut State,
)
    -> three_d::FrameOutput
{
    update_gui_and_viewport(frame_input, gui, state);

    rust_dsymbols::display::controls::orbit_control_update_camera(
        &mut state.camera, &mut frame_input.events, 1.0, 1000.0
    );

    let [r, g, b, a] = state.options.background_color.to_normalized_gamma_f32();
    let clear_state = three_d::ClearState::color_and_depth(r, g, b, a, 1.0);
    let screen = frame_input.screen();

    state.update_caches();

    if let Some(parts) = state.parts_cache_key()
        .and_then(|k| state.parts_cache.get(&k))
    {
        match parts {
            Ok(parts) => {
                state.message = String::from("Ok!");

                let options = &state.options;
                let models = build_models(context, parts, options);
                screen.clear(clear_state);

                render_models(context, &screen, options, &state.camera, models);
            },
            Err(msg) => {
                state.message = msg.clone();
                screen.clear(clear_state);
            },
        }
    } else {
        state.message = String::from("nothing to render");
        screen.clear(clear_state);
    }

    screen.write(|| gui.render()).unwrap();

    // Ensures a valid return value.
    three_d::FrameOutput::default()
}


fn render_models(
    context: &mut three_d::Context,
    screen: &three_d::RenderTarget<'_>,
    options: &Options,
    camera: &three_d::Camera,
    models: Vec<three_d::Gm<three_d::InstancedMesh, three_d::PhysicalMaterial>>
) {
    let sun_dir = (
        camera.view().invert().unwrap() *
        options.sun_direction.extend(0.0)
    ).truncate();

    let white = three_d::Srgba::WHITE;
    let mut sun = three_d::DirectionalLight::new(context, 2.0, white, sun_dir);
    let ambient = three_d::AmbientLight::new(context, 0.1, white);

    if options.sun_casts_shadows {
        sun.generate_shadow_map(9192, &models);
    }
    screen.render(camera, models, &[&sun, &ambient]);
}


fn update_gui_and_viewport(
    frame_input: &mut three_d::FrameInput,
    gui: &mut three_d::GUI,
    state: &mut State)
{
    let mut offset_wd = 0.0;
    let mut offset_ht = 0.0;

    gui.update(
        &mut frame_input.events,
        frame_input.accumulated_time,
        frame_input.viewport,
        frame_input.device_pixel_ratio,
        |gui_context| {
            (offset_wd, offset_ht) = gui_callback(state, gui_context);
        },
    );

    let wd = offset_wd * frame_input.device_pixel_ratio;
    let ht = offset_ht * frame_input.device_pixel_ratio;

    let viewport = three_d::Viewport {
        x: wd as i32,
        y: ht as i32,
        width: frame_input.viewport.width - wd as u32,
        height: frame_input.viewport.height - ht as u32,
    };
    state.camera.set_viewport(viewport);
}


fn gui_callback(state: &mut State, gui_context: &three_d::egui::Context)
    -> (f32, f32)
{
    three_d::egui::SidePanel::left("side_panel").show(gui_context, |ui| {
        ui.heading("Tiling");
        ui.add_space(12.0);

        ui.label(format!("Collection '{}'", state.catalog.collection_name));

        ui_file_loader(ui, gui_context, state);
        ui.add_space(12.0);

        ui_navigation_buttons(ui, state);
        ui.add_space(24.0);

        ui.heading("Appearance");
        ui.add_space(12.0);

        let options = &mut state.options;

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

    let width = gui_context.used_rect().width();

    let response = three_d::egui::TopBottomPanel::top("status line").show(
        gui_context,
        |ui| { ui.label(&state.message[..]); }
    );
    let height = response.response.rect.height();

    (width, height)
}


fn ui_file_loader(
    ui: &mut three_d::egui::Ui,
    gui_context: &three_d::egui::Context,
    state: &mut State,
) {
    if (ui.button("Load...")).clicked() {
        let filter = Box::new({
            let ext = Some(std::ffi::OsStr::new("ds"));
            move |path: &Path| -> bool { path.extension() == ext }
        });
        let mut dialog = FileDialog::open_file(state.opened_file.clone())
            .show_files_filter(filter);
        dialog.open();
        state.open_file_dialog = Some(dialog);
    }

    if let Some(dialog) = &mut state.open_file_dialog {
        if dialog.show(gui_context).selected() {
            if let Some(file) = dialog.path() {
                let path = file.to_path_buf();
                state.opened_file = Some(path.clone());
                match file_contents(&path) {
                    Ok(content) => {
                        let title = path.components().last().unwrap()
                            .as_os_str().to_str().unwrap();
                        let collection = content.lines().map(String::from);
                        state.catalog.add_collection(title, collection);
                    },
                    Err(err) => {
                        state.message = err.to_string()
                    },
                }
            }
        }
    }
}


fn ui_navigation_buttons(
    ui: &mut three_d::egui::Ui,
    state: &mut State,
) {
    let n = state.catalog.collections[&state.catalog.collection_name].len();
    let k = state.catalog.index_in_collection;

    ui.label(format!("Index {} of {}", k + 1, n));

    ui.horizontal(|ui| {
        if k == 0 {
            ui.add_enabled(false, three_d::egui::Button::new("First"));
            ui.add_enabled(false, three_d::egui::Button::new("Prev"));
        } else {
            if ui.button("First").clicked() {
                state.catalog.index_in_collection = 0;
            }
            if ui.button("Prev").clicked() {
                state.catalog.index_in_collection -= 1;
            }
        }
        if k + 1 >= n {
            ui.add_enabled(false, three_d::egui::Button::new("Next"));
            ui.add_enabled(false, three_d::egui::Button::new("Last"));
        } else {
            if ui.button("Next").clicked() {
                state.catalog.index_in_collection += 1;
            }
            if ui.button("Last").clicked() {
                state.catalog.index_in_collection = n.max(1) - 1;
            }
        }
    });
}


fn file_contents(path: &PathBuf) -> std::io::Result<String> {
    let mut file = File::open(path)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    Ok(contents)
}


fn build_models(
    context: &three_d::Context,
    parts: &Vec<(three_d::CpuMesh, ItemType)>,
    options: &Options
)
    -> Vec<three_d::Gm<three_d::InstancedMesh, three_d::PhysicalMaterial>>
{
    let instances = three_d::Instances {
        transformations: vec![
            Mat4::from_scale(1.0),
        ],
        ..Default::default()
    };

    parts.iter().map(|(part_mesh, item_type)| {
        let color = match item_type {
            ItemType::Vertex => options.vertex_color,
            ItemType::Edge => options.edge_color,
            ItemType::Face => options.face_color,
        };

        three_d::Gm::new(
            three_d::InstancedMesh::new(context, &instances, part_mesh),
            three_d::PhysicalMaterial {
                albedo: color.to_normalized_gamma_f32().into(),
                metallic: 0.0,
                roughness: 0.5,
                ..Default::default()
            }
        )
    }).collect()
}


fn build_parts(til: &Tiling, options: &Options)
    -> Result<Vec<(three_d::CpuMesh, ItemType)>, String>
{
    let mut parts = vec!();

    for tile_mesh in tiles(til)? {
        for (part_mesh, item_type) in decompose_mesh(
            &scaled_mesh(&tile_mesh, options.tile_scale),
            options.edge_radius
        ) {
            parts.push((part_mesh.to_cpu_mesh(), item_type))
        }
    }

    Ok(parts)
}


fn tiles(til: &Tiling)
    -> Result<Vec<Mesh<Point3<f64>>>, String>
{
    let ds = til.ds().as_ref()?;
    let cov = til.cov();
    let cov = cov.as_ref().as_ref()?;
    let skel = til.skel();
    let skel = skel.as_ref().as_ref()?;

    let basis = invariant_basis(
        &gram_matrix(ds, cov, skel).ok_or("error computing the Gram matrix")?
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

        Mesh::from_oriented_faces(vs, faces.clone())
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
