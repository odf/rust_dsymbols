use nom::branch::alt;
use nom::character::complete::{char, digit1, multispace0, multispace1};
use nom::multi::separated_list1;
use nom::sequence::{separated_pair, tuple};
use nom::IResult;
use nom::combinator::{map, map_opt};


#[derive(Debug, Eq, PartialEq)]
pub struct DSymSpec {
    pub set_count: usize,
    pub sym_count: usize,
    pub size: usize,
    pub dim: usize,
    pub op_spec: Vec<Vec<usize>>,
    pub m_spec: Vec<Vec<usize>>
}


pub fn parse_dsymbol(input: &str) -> Result<(&str, DSymSpec), String>
{
    dsymbol(input).map_err(|e| e.to_string())
}


fn dsymbol(input: &str) -> IResult<&str, DSymSpec> {
    map(
        tuple((
            tuple((ws0, char('<'), ws0)),
            counts,
            tuple((ws0, char(':'), ws0)),
            extents,
            tuple((ws0, char(':'), ws0)),
            int_lists,
            tuple((ws0, char(':'), ws0)),
            int_lists,
            tuple((ws0, char('>'), ws0)),
        )),
        |(_, (set_count, sym_count), _, (size, dim), _, op_spec, _, m_spec, _)|
        {
            DSymSpec { set_count, sym_count, size, dim, op_spec, m_spec }
        }
    )(input)
}


fn counts(input: &str) -> IResult<&str, (usize, usize)> {
    separated_pair(integer, char('.'), integer)(input)
}


fn extents(input: &str) -> IResult<&str, (usize, usize)> {
    alt((
        separated_pair(integer, ws1, integer),
        map(integer, |n| (n, 2))
    ))(input)
}


fn int_lists(input: &str) -> IResult<&str, Vec<Vec<usize>>> {
    separated_list1(tuple((ws0, char(','), ws0)), int_list)(input)
}


fn int_list(input: &str) -> IResult<&str, Vec<usize>> {
    separated_list1(ws1, integer)(input)
}


fn integer(input: &str) -> IResult<&str, usize> {
    map_opt(digit1, |digits: &str| digits.parse().ok())(input)
}


fn ws0(input: &str) -> IResult<&str, &str> {
    multispace0(input)
}


fn ws1(input: &str) -> IResult<&str, &str> {
    multispace1(input)
}


#[test]
fn test_parse_tiny_dsymbol() {
    assert_eq!(
        dsymbol("<1.1:1:1,1,1:3,4>"),
        Ok(("", DSymSpec {
            set_count: 1,
            sym_count: 1,
            size: 1,
            dim: 2,
            op_spec: vec![vec![1], vec![1], vec![1]],
            m_spec: vec![vec![3], vec![4]],
        }))
    );
}


#[test]
fn test_parse_larger_dsymbol() {
    assert_eq!(
        dsymbol("<10.8:2 3:1 2,1 2,1 2,2:3 3,3 4,4>"),
        Ok(("", DSymSpec {
            set_count: 10,
            sym_count: 8,
            size: 2,
            dim: 3,
            op_spec: vec![vec![1, 2], vec![1, 2], vec![1, 2], vec![2]],
            m_spec: vec![vec![3, 3], vec![3, 4], vec![4]],
        }))
    );
}


#[test]
fn test_parse_incomplete_dsymbol() {
    assert_eq!(
        dsymbol("<1.1:1:1,1,0:3,0>"),
        Ok(("", DSymSpec {
            set_count: 1,
            sym_count: 1,
            size: 1,
            dim: 2,
            op_spec: vec![vec![1], vec![1], vec![0]],
            m_spec: vec![vec![3], vec![0]],
        }))
    );
}


#[test]
fn test_parse_with_extra_spaces() {
    assert_eq!(
        dsymbol(" < 10.8: 2 3:1 2 , 1 2,  1 2  ,2 :3 3 , 3   4,4 >  "),
        Ok(("", DSymSpec {
            set_count: 10,
            sym_count: 8,
            size: 2,
            dim: 3,
            op_spec: vec![vec![1, 2], vec![1, 2], vec![1, 2], vec![2]],
            m_spec: vec![vec![3, 3], vec![3, 4], vec![4]],
        }))
    );
}


#[test]
fn test_parse_multiline() {
    assert_eq!(
        dsymbol(
            "<10.8:2 3:
            1 2,1 2,1 2,2
            :3 3,3 4,4>"
        ),
        Ok(("", DSymSpec {
            set_count: 10,
            sym_count: 8,
            size: 2,
            dim: 3,
            op_spec: vec![vec![1, 2], vec![1, 2], vec![1, 2], vec![2]],
            m_spec: vec![vec![3, 3], vec![3, 4], vec![4]],
        }))
    );
}


#[test]
fn test_parse_malformed() {
    assert!(dsymbol("").is_err());
    assert!(dsymbol("<>").is_err());
    assert!(dsymbol("1.1:1:1,1,1:3,4").is_err());
    assert!(dsymbol("<<1.1:1:1,1,1:3,4>>").is_err());
    assert!(dsymbol("<1. 1:1:1,1,1:3,4>").is_err());
    assert!(dsymbol("<1.1:1:1,1,1:3,4:>").is_err());
    assert!(dsymbol("<1.1::1,1,1:3,4>").is_err());
    assert!(dsymbol("<1.1:1:1,,1:3,4>").is_err());
}
