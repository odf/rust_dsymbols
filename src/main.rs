use rust_dsymbols::backtrack_demo::IntPartitions;


fn main() {
    for xs in IntPartitions::new(10) {
        println!("{:?}", xs);
    }
}
