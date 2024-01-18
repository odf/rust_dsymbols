use nom::branch::alt;
use nom::character::complete::{char, digit1, space1, space0};
use nom::multi::separated_list1;
use nom::sequence::{separated_pair, tuple};
use nom::{IResult, Finish};
use nom::combinator::{map, map_opt};


pub struct DSymSpec {
    pub set_count: u32,
    pub sym_count: u32,
    pub size: u32,
    pub dim: u32,
    pub op_spec: Vec<Vec<u32>>,
    pub v_spec: Vec<Vec<u32>>
}


pub fn parse_dsymbol(input: &str) -> Result<(&str, DSymSpec), String>
{
    dsymbol(input).finish().map_err(|e| e.to_string())
}


fn dsymbol(input: &str) -> IResult<&str, DSymSpec> {
    map(
        tuple((
            tuple((space0, char('<'), space0)),
            counts,
            tuple((space0, char(':'), space0)),
            extents,
            tuple((space0, char(':'), space0)),
            int_lists,
            tuple((space0, char(':'), space0)),
            int_lists,
            tuple((space0, char('>'), space0)),
        )),
        |(_, (set_count, sym_count), _, (size, dim), _, op_spec, _, v_spec, _)|
        {
            DSymSpec { set_count, sym_count, size, dim, op_spec, v_spec }
        }
    )(input)
}


fn counts(input: &str) -> IResult<&str, (u32, u32)> {
    separated_pair(integer, char('.'), integer)(input)
}


fn extents(input: &str) -> IResult<&str, (u32, u32)> {
    alt((
        separated_pair(integer, space1, integer),
        map(integer, |n| (n, 2))
    ))(input)
}


fn int_lists(input: &str) -> IResult<&str, Vec<Vec<u32>>> {
    separated_list1(tuple((space0, char(','), space0)), int_list)(input)
}


fn int_list(input: &str) -> IResult<&str, Vec<u32>> {
    separated_list1(space1, integer)(input)
}


fn integer(input: &str) -> IResult<&str, u32> {
    map_opt(digit1, map_integer)(input)
}


fn map_integer(digits: &str) -> Option<u32> {
    if digits.len() <= 9 {
        Some(digits.chars().fold(0, |n, c| n * 10 + c.to_digit(10).unwrap()))
    } else {
        None
    }
}
